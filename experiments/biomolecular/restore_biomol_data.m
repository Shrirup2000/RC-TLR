function restore_biomol_data(dataRoot,names)
%RESTORE_BIOMOL_DATA Restore the exact original MAT files from binary parts.
% Existing files are never overwritten; distributed data are SHA-256 checked.
% This routine uses the Java runtime bundled with standard desktop MATLAB.
    here=fileparts(mfilename('fullpath'));
    if nargin<1 || isempty(dataRoot), dataRoot=fullfile(here,'data'); end
    manifest=fullfile(dataRoot,'manifest.csv');
    if ~isfile(manifest), return; end % External, user-provided data folder.
    T=readtable(manifest,'TextType','string');
    if nargin<2, names=T.File; end
    names=string(names);
    for i=1:numel(names)
        row=T(T.File==names(i),:);
        assert(height(row)==1,'Dataset missing from manifest: %s',names(i));
        found=dir(fullfile(dataRoot,'**',char(row.File)));
        if ~isempty(found)
            assert(numel(found)==1,'Multiple copies of %s',row.File);
            path=fullfile(found.folder,found.name);
            check_file(path,row);
            continue;
        end
        assert(row.Parts>0,'Missing bundled data file: %s',row.File);
        target=fullfile(dataRoot,char(row.File)); tmp=[tempname(dataRoot),'.tmp'];
        cleanup=onCleanup(@()remove_temp(tmp)); %#ok<NASGU>
        join_parts(tmp,dataRoot,row);
        check_file(tmp,row);
        assert(~isfile(target),'Refusing to overwrite %s',target);
        [ok,msg]=movefile(tmp,target);
        assert(ok,'Could not restore %s: %s',target,msg);
        fprintf('Restored and verified %s\n',row.File);
        clear cleanup
    end
end

function join_parts(path,root,row)
    out=fopen(path,'wb'); assert(out~=-1,'Cannot create %s',path);
    closeOut=onCleanup(@()fclose(out)); %#ok<NASGU>
    for k=1:row.Parts
        name=sprintf('%s.part%03d',row.File,k);
        part=fullfile(root,'parts',name);
        in=fopen(part,'rb'); assert(in~=-1,'Missing binary dataset part: %s',part);
        closeIn=onCleanup(@()fclose(in));
        bytes=fread(in,inf,'*uint8');
        assert(fwrite(out,bytes,'uint8')==numel(bytes),'Incomplete write of %s',name);
        clear closeIn
    end
end

function check_file(path,row)
    info=dir(path);
    assert(info.bytes==row.Bytes,'Byte count mismatch for %s',row.File);
    md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
    in=fopen(path,'rb'); assert(in~=-1,'Cannot read %s',path);
    closeIn=onCleanup(@()fclose(in)); %#ok<NASGU>
    while true
        bytes=fread(in,1048576,'*uint8');
        if isempty(bytes), break; end
        md.update(typecast(bytes(:),'int8'));
    end
    digest=typecast(md.digest(),'uint8');
    hex=lower(reshape(dec2hex(digest,2).',1,[]));
    assert(strcmp(hex,row.SHA256),'SHA-256 mismatch for %s',row.File);
end

function remove_temp(path)
    if isfile(path), delete(path); end
end
