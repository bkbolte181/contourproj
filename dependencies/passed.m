function [was_passed] = passed(parameter, arguments)
%PASSED Check if an argument was passed to vargin
% Parameters
%   argument: varagin
%   parameter: parameter which was passed

v = deal(arguments);

was_passed = 0;

for i=1:numel(v)
    a = v{i};
    if a == parameter
        was_passed = 1;
        break;
    end
end

end

