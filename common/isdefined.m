function r = isdefined(var_s, vars)
%% Cheacks whether a variable is initialized

if any(strcmp(var_s, vars)) 
    r = ~isempty(evalin('caller',var_s));
else
    r = false;
end