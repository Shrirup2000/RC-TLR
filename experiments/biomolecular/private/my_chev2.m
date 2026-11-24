function X = my_chev2(n)
%MY_CHEV2 First-kind Chebyshev nodes on [-1,1], in descending order.
    validateattributes(n, {'numeric'}, {'scalar','integer','positive','finite'});
    X = cos((2*(1:double(n))-1)*(pi/(2*double(n))));
end
