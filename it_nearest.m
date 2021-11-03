function c =  it_nearest(vec,num)

% returns the index of the element in vec which is closest to num

a = find(vec>=num);
if isempty(a)
    a = length(vec);
end
b = find(vec<=num);
if isempty(b)
    b = 1;
end
c = [];
while isempty(c)
    c = intersect(a,b);
    b(end) = b(end) +1;
end