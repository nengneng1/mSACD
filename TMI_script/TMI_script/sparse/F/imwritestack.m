function imwritestack(stack, filename)

t = Tiff(filename, 'w');

tagstruct.ImageLength = size(stack, 1);
tagstruct.ImageWidth = size(stack, 2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

for k = 1:size(stack, 3)
    t.setTag(tagstruct)
    t.write(uint16(stack(:, :, k)));
    t.writeDirectory();
    if mod(k,100)==0
        k;
    end
end

t.close();
end