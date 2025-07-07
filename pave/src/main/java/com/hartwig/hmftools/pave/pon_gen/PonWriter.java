package com.hartwig.hmftools.pave.pon_gen;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.pon_gen.PonConfig.GERMLINE_CLINVAR_MAX_REPEAT;
import static com.hartwig.hmftools.pave.pon_gen.PonConfig.GERMLINE_CLINVAR_MIN_SAMPLES;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class PonWriter
{
    private final PonConfig mConfig;

    private final BufferedWriter mWriter;
    private final List<VariantPonData> mManualEntries;

    private final List<ChrBaseRegion> mRemainingRegions;
    private final Map<ChrBaseRegion,List<VariantPonData>> mCachedRegionVariants;

    public PonWriter(final PonConfig config, final List<ChrBaseRegion> regions, final List<VariantPonData> manualEntries)
    {
        mConfig = config;
        mManualEntries = manualEntries;

        mWriter = initialiseWriter();

        mRemainingRegions = Lists.newArrayList(regions);
        mCachedRegionVariants = Maps.newHashMap();
    }

    public synchronized void onVariantsComplete(final ChrBaseRegion region, final List<VariantPonData> variants)
    {
        // now include any manual entries from this region
        mManualEntries.stream().filter(x -> region.containsPosition(x.Chromosome, x.Position)).forEach(x -> variants.add(x));

        Collections.sort(variants, new VariantPonData.VariantSorter());

        mCachedRegionVariants.put(region, variants);
        writeCompleteRegions();
    }

    public void close()
    {
        if(!mCachedRegionVariants.isEmpty())
        {
            PV_LOGGER.error("closing PON writer with {} cached variants", mCachedRegionVariants.size());

            for(ChrBaseRegion region : mRemainingRegions)
            {
                List<VariantPonData> variants = mCachedRegionVariants.remove(region);

                if(variants != null)
                {
                    writeVariants(variants);
                }
            }
        }

        closeBufferedWriter(mWriter);
    }

    private void writeCompleteRegions()
    {
        while(!mRemainingRegions.isEmpty())
        {
            ChrBaseRegion region = mRemainingRegions.get(0);

            List<VariantPonData> variants = mCachedRegionVariants.remove(region);

            if(variants == null)
                return;

            mRemainingRegions.remove(0);

            writeVariants(variants);
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            // PV_LOGGER.info("writing {} variants to PON file({})", variants.size(), fileName);

            BufferedWriter writer = createBufferedWriter(mConfig.OutputFilename);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME);
            sj.add(FLD_POSITION);
            sj.add(FLD_REF);
            sj.add(FLD_ALT);

            sj.add("SampleCount");
            sj.add("MaxSampleReadCount");
            sj.add("TotalReadCount");

            if(!mConfig.WriteFinal)
            {
                sj.add("SomaticHotspot");
                sj.add("GermlineHotspot");
                sj.add("ClinvarPathogenic");
                sj.add("InCodingRegion");
                sj.add("RepeatCount");
            }

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise output file: {}", e.toString());
            System.exit(1);
            return null;
        }
    }

    private void writeVariants(final List<VariantPonData> variants)
    {
        try
        {
            for(VariantPonData variant : variants)
            {
                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(variant.Chromosome);
                sj.add(String.valueOf(variant.Position));
                sj.add(variant.Ref);
                sj.add(variant.Alt);
                sj.add(String.valueOf(variant.sampleCount()));
                sj.add(String.valueOf(variant.maxSampleReadCount()));
                sj.add(String.valueOf(variant.totalReadCount()));

                if(mConfig.WriteFinal)
                {
                    if(filterOutVariant(variant))
                        continue;
                }
                else
                {
                    sj.add(String.valueOf(variant.isSomaticHotspot()));
                    sj.add(String.valueOf(variant.isGermlineHotspot()));
                    sj.add(String.valueOf(variant.isClinvarPathogenic()));
                    sj.add(String.valueOf(variant.inCodingRegion()));
                    sj.add(String.valueOf(variant.repeatCount()));
                }

                mWriter.write(sj.toString());
                mWriter.newLine();
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise output file: {}", e.toString());
            System.exit(1);
        }
    }

    private boolean filterOutVariant(final VariantPonData variant)
    {
        // exckude a variant from the PON if either of the following dot points are satisfied:
        // - somatic hotspot
        // - (germline hotspot OR clinvar pathogenic) AND has sampleCount < 10 AND (has repeatCount < 4 or is NOT an indel)
        if(variant.isSomaticHotspot())
            return true;

        if(variant.sampleCount() < GERMLINE_CLINVAR_MIN_SAMPLES
        && (variant.isGermlineHotspot() || variant.clinvarPathogenicity().isPathogenic())
        && (!variant.isIndel() || variant.repeatCount() < GERMLINE_CLINVAR_MAX_REPEAT))
        {
            return true;
        }

        return false;

    }
}
