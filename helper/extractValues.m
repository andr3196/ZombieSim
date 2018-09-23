function str = extractValues(filename)
    text = fileread(filename);
    b = regexp(text, '%%% PARAM BEGIN %%%', 'once');
    e = regexp(text, '%%% PARAM END %%%', 'once');
    str = text(b + 19: e-1);
