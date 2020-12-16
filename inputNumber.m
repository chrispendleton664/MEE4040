function out = inputNumber(prompt)
out = NaN;
tot = 3; % how many attempts to make
while isnan(out) && tot>0
    out = str2double(input('Write your age :', 's'));
    tot = tot-1;
end
end