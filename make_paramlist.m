function ret = make_paramlist(paramvec)
    % 生成参数列表
    dc = paramvec(1);
    df = paramvec(2);
    rr = paramvec(3);
    grf = paramvec(4);
    dmi = paramvec(5);
    lis = paramvec(6);
    red = paramvec(7);
    ret = struct( ...
    'drive_conversion', dc, ...
    'drive_fitness', df, ...
    'release_rate', rr, ...
    'germline_resistance_forming', grf, ...
    'dd_mothertal_inheritance', dmi,...
    'log_immunity_speed', lis,...
    'reducehtm',red ...
        );
end