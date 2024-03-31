        function frame2gif(fr, filename, varargin, kwargs, im_args)
            % FRAME2GIF - Write an array of frames to a GIF file
            %
            % FRAME2GIF(fr, filename) writes the series of frames fr
            % to the file filename.
            %
            % FRAME2GIF(fr) uses 'tmp.gif', as the default filename
            %
            % FRAME2GIF(...,'map', N) creates a colorspace with N values.
            % The defaul is 2^8.
            %
            % FRAME2GIF(..., Name, Value, ...) forwards Name/Value pairs to
            % imwrite.
            %
            % Example:
            % % create some data
            % sz = [151,101,100];
            % im = complex(rand(sz), rand(sz)) - (0.5 + 0.5i);
            % im = rad2deg(angle(im)); % display phase in degrees
            % 
            % % display with imagesc
            % hf = figure();
            % him = imagesc(im(:,:,1));
            % colormap hsv;
            % colorbar;
            % title('Random phase');
            % 
            % % collect the frames
            % fps = 10; % frame per second
            % getframe(hf); % work-around: throw away 1 frame to initialize
            % for i = 1:size(im,3)
            %     him.CData(:) = im(:,:,i);
            %     drawnow; % pause(1/fps);
            %     fr(i) = getframe(hf);
            % end
            %
            % % Write frames into a GIF
            % frame2gif(fr, 'noise', 'LoopCount', Inf, 'DelayTime', 1/fps);
            %
            % % Customize the GIF with valid arguments for imwrite
            % frame2gif(fr, 'WriteMode', 'append', 'DisposalMethod', 'restorePrevious')
            % 
            % See also ANIMATE SCAN/GIF CHANNELDATA/GIF IMWRITE

            arguments
                fr struct
                filename (1,1) string = 'tmp.gif'
            end
            arguments(Repeating)
                varargin % unrecognized inputs (forwarded to imwrite)
            end
            arguments
                kwargs.dither (1,1) logical = false
                kwargs.map (:,:) double = 2^8
                im_args.WriteMode (1,1) string {mustBeMember(im_args.WriteMode, ["overwrite","append"])}
                im_args.LoopCount (1,1) double {mustBePositive} = Inf
                im_args.DelayTime (1,1) double {mustBePositive} = 1/15
                im_args.Comment (1,:) string
                im_args.DisposalMethod (1,1) string {mustBeMember(im_args.DisposalMethod, ["leaveInPlace","restoreBG","restorePrevious","doNotSpecify"])}
                im_args.TransparentColor (1,1) double {mustBeInteger}
                im_args.BackgroundColor (1,1) double {mustBeInteger}
                im_args.ScreenSize (1,2) double 
                im_args.Location (1,2) double 
            end

            % parse inputs
            for i = 1:2:numel(varargin)
                switch varargin{1}
                    case {'map', 'dither'}
                        kwargs.(varargin{1}) = varargin{i+1};
                    otherwise
                        im_args.(varargin{1}) = varargin{i+1};
                end
            end

            % enforce filename ending
            if ~endsWith(filename, '.gif','IgnoreCase',true)
                filename = filename + ".gif";
            end

            if kwargs.dither, opt_dither = 'dither'; else, opt_dither = 'nodither'; end

            % get color space for the image
            [~, map] = rgb2ind(fr(1).cdata, kwargs.map, opt_dither);

            % get all image data
            im = arrayfun(@(fr) {rgb2ind(fr.cdata,map,opt_dither)}, fr);
            im = cat(4, im{:});

            % forward valid options to imwrite 
            nvkwargs = struct2nvpair(im_args);
            imwrite(im, map, filename, nvkwargs{:});
        end