function [sonout,Fs] = macg(targetname, bricknames, outfname)
% Matlab function used to create the audio in the Youtube video
% "Audio shenanigans with deconvolution," on the GhostBass channel.

% "target" should be the filename of a sound to compose.
% "bricks" should be a cell array of filenames of sounds to use in
% composing the target.

% Example call:
% [so,Fs] = macg('ranggyu.mp3', {'dog1.wav','dog2.wav','dog3.wav','dog4.wav','dog5.wav','dog6.wav'}, 'roll.wav');

% Read in the audio of the original song:
[sont,Fs] = audioread(targetname);

% We're just going to output mono, so average the two channels if the input is stereo:
if size(sont,2)>1
    sont = sum(sont,2);
end

oN = length(sont);  %  original number of samples
N = 2^floor(log(oN)/log(2));  %  cut down to a power of 2 for easier fft
sont = sont(oN+1-N:oN);
fftt = fft(sont);

% Now load in the dog bark samples:
B = length(bricknames);
brix = zeros(N,B);
maxcobblelength = 0;
for b = 1:B
    [brin,fsin] = audioread(bricknames{b});
    if size(brin,2)>1
        brin = sum(brin,2);  %  will scale result at the end
    end
    % Fix any discrepancy in sample rate (though actually they were all 44100 Hz):
    if fsin ~= Fs
        brin = resample(timeseries(brin,(1:length(brin))/fsin), 1/Fs:1/Fs:length(brin)/fsin);
    end
    maxcobblelength = max([maxcobblelength,length(brin)]);  %  length of longest bark
    brix(1:length(brin),b) = brin;
end
fftb = fft(brix);  %  columnwise

% In principle, we only need one brick, but using a few may sound more
% natural, and as a side benefit, it allows us to avoid very high ratios
% between the target's fft and the bricks' ffts. So here we split the
% original signal's fft into B parts, load-balanced across the bricks.
fk = (fftt ./ sum(abs(fftb),2)) .* abs(fftb);
% So now each column of fk will be built from the corresponding brick.

Int = ifft(fk ./ fftb);
% Now each column, convolved with its brick (the corresponding column in
% brix), should reconstitute the corresponding part of the target. (So to
% reconstruct the original mono track exactly, we could just use conv to
% combine each column of Int with the corresponding brick, and add them
% up.) We will do a "pointilised" imperfect convolution, with resolution
% gradually increasing over the target's duration.

% Initially, I had
  % intensity = exp(linspace(log(2/(B*Fs)), log(1e8), N))';
% but actually that was too fast toward the end. So jockey the Poisson intensity for a more dramatic effect:
intensity = exp(polyval(polyfit([0,0.05,0.25,0.5,0.75,1],[log(1/(6*44100)),-11.5,-2.5,2.5,5.5,log(1e7)],5),linspace(0,1,N)))';

sonout = zeros(N+maxcobblelength-1,1);
for b=1:B
    convand = Int(:,b) .* poissrnd(intensity) .* intensity.^(-2/3);  %  it seems
    convand(end) = Int(end,b);  %  otherwise =0*inf=NaN
    sonout = sonout + cconv(convand,brix(:,b),N+maxcobblelength-1);
end

% At this point we've done all the real work; the rest is just "mastering."

sonout = sonout(1:N);  %  remove extra silence at the end

% Now use kernel smoothing to approximately level the volume across the track:
kerhw = floor(Fs/4); %  force kernel width to be odd
epker = (1:2*kerhw+1) .* (2*kerhw+1:-1:1); epker = epker';
volenv = conv(abs(sonout), epker);  %  volume envelope
volenv = volenv + 1e-6*max(volenv);  %  to avoid 0/0 noise (especially at the beginning)
sonout = sonout ./ volenv(kerhw+1:end-kerhw);

% Tidy the beginning and end:
sonout(1:kerhw) = sonout(1:kerhw) .* linspace(0,1,kerhw)';  %  quick fade in to hide edge effects of kernel smoothing
sonout(end-15*Fs+1:end) = sonout(end-15*Fs+1:end) .* linspace(1,0,15*Fs)'.^2;  %  fade out at the end (similar shape to original song)

% Normalise the overall volume and output the audio:
sonout = sonout / max(abs(sonout));
audiowrite(outfname, sonout, Fs);
