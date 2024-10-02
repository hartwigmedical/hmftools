package com.hartwig.hmftools.chord.variant;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.io.File;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.variant.variantcontext.VariantContext;

public class VcfFile
{
    public final String mPath;
    private final boolean mIncludeNonPass;

    public VcfFile(String path, boolean includeNonPass)
    {
        mPath = path;
        mIncludeNonPass = includeNonPass;
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
        if(variants.size()==0)
        {
            CHORD_LOGGER.warn("No {}variants found in vcf file: {}", filterType, mPath);
        }
        else
        {
            CHORD_LOGGER.info("Loaded {} {}variants from vcf file: {}", variants.size(), filterType, mPath);
        }

        return variants;
    }
}
