function [peaks, nNeurons, nFrames] = fn_loadPeaksData(dataFile, chunkSize)
% FN_LOADPEAKSDATA  Chunked load of a large logical 'peaks' matrix from disk.
%
%   [peaks, nNeurons, nFrames] = FN_LOADPEAKSDATA(dataFile, chunkSize)
%
%   Loads the variable 'peaks' (nNeurons x nFrames) from a .mat file using
%   matfile() partial-loading, reading CHUNKSIZE neurons (rows) at a time
%   to keep peak memory usage low for very large recordings. The matrix on
%   disk is assumed to be numeric (double/logical); it is cast to logical
%   in memory, which is typically an order of magnitude smaller.
%
% INPUTS
%   dataFile  - path to a .mat file containing a variable 'peaks'
%   chunkSize - (optional) number of neurons to load per chunk (default 50)
%
% OUTPUTS
%   peaks     - nNeurons x nFrames logical spike/event matrix
%   nNeurons  - number of neurons (rows)
%   nFrames   - number of frames (columns)

if nargin < 2 || isempty(chunkSize)
    chunkSize = 50;
end

fprintf('  Loading %s ...\n', dataFile);
mf   = matfile(dataFile);
info = whos(mf, 'peaks');
info = info(1);
sz   = info.size;
nNeurons = sz(1);
nFrames  = sz(2);

fprintf('  File: %d neurons x %d frames\n', nNeurons, nFrames);
fprintf('  On disk: %.1f GB (double)  ->  logical in memory: %.0f MB\n', ...
    info.bytes/1e9, nNeurons*nFrames/1e6);

peaks    = false(nNeurons, nFrames);
n_chunks = ceil(nNeurons / chunkSize);
fprintf('  Loading %d chunks of %d neurons...\n', n_chunks, chunkSize);

for c = 1:n_chunks
    r0 = (c-1)*chunkSize + 1;
    r1 = min(c*chunkSize, nNeurons);
    peaks(r0:r1,:) = logical(mf.peaks(r0:r1,:));
    if mod(c,5)==0 || c==n_chunks
        fprintf('    chunk %d/%d  (neurons %d-%d)\n', c, n_chunks, r0, r1);
    end
end

end
