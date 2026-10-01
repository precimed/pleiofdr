function save_figure(handle, format, filename)
    if isempty(handle), return; end;
    if length(handle) > 1, handle=handle(1); end;
    screen_size = get(0, 'ScreenSize');
    origSize = get(handle, 'Position'); % grab original on screen size
    set(handle, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
    %set(handle,'PaperPositionMode','auto') %set paper pos for printing
    set(handle,'PaperPositionMode','manual') %set paper pos for printing

    if (strcmp(format, 'png') && ~is_octave())
        if isobject(handle), handle_Number = handle.Number; else handle_Number = handle; end;
        print(filename,'-dpng','-r0',sprintf('-f%i', handle_Number));
    else
        outfn = sprintf('%s.%s', filename, format);
        fprintf('writing %s\n', outfn);
        saveas(handle, outfn, format);
    end

    set(handle, 'Position', origSize) %set back to original dimensions
end
