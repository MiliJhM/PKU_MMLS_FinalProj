function make_gif(filename, index)
    f = getframe(gcf);
    imind = frame2im(f);
    [imind, cm] = rgb2ind(imind,256);
    if index == 1
        imwrite(imind, cm, filename, 'gif', 'LoopCount', inf, 'DelayTime', 0.0001);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode','append','DelayTime', 0.0001);
    end
end