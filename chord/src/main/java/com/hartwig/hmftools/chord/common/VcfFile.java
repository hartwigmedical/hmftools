package com.hartwig.hmftools.chord.common;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.io.File;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.variant.variantcontext.VariantContext;

public class VcfFile implements LoggingOptions
{
    public final String mPath;
    public final boolean mIncludeNonPass;
    public String mLogPrefix = "";

    public VcfFile(String path, boolean includeNonPass)
    {
        mPath = path;
        mIncludeNonPass = includeNonPass;
    }

    @Override
    public VcfFile logPrefix(String logPrefix)
    {
        mLogPrefix = logPrefix;
        return this;
    }

    public List<VariantContext> loadVariants() throws NoSuchFileException
    {
        if(!new File(mPath).isFile())
            throw new NoSuchFileException(mPath);

        VcfFileReader vcfFileReader = new VcfFileReader(mPath);

        List<VariantContext> variants = new ArrayList<>();

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            boolean isPassVariant = variantContext.getFilters().isEmpty();

            if(!mIncludeNonPass && !isPassVariant)
                continue;

            variants.add(variantContext);
        }

        String filterType = (mIncludeNonPass) ? "" : "PASS ";
        CHORD_LOGGER.debug("{}Loaded {} {}variants from: {}", mLogPrefix, variants.size(), filterType, mPath);

        return variants;
    }
}
