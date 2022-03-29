package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class Blacklistings
{
    private final List<ChrBaseRegion> mBedRegions;
    private final List<VcfEntry> mVcfEntries;
    private boolean mHasValidData;

    private static final String BLACKLIST_BED = "blacklist_bed";
    private static final String BLACKLIST_VCF = "blacklist_vcf";

    private static final String BLACKLIST_BED_FLAG = "BLACKLIST_BED";
    private static final String BLACKLIST_VCF_FLAG = "BLACKLIST_VCF";

    public Blacklistings(final CommandLine cmd)
    {
        mBedRegions = Lists.newArrayList();
        mVcfEntries = Lists.newArrayList();
        mHasValidData = true;

        if(cmd.hasOption(BLACKLIST_BED))
        {
            loadBedEntries(cmd.getOptionValue(BLACKLIST_BED));
        }

        if(cmd.hasOption(BLACKLIST_VCF))
        {
            loadVcfEntries(cmd.getOptionValue(BLACKLIST_VCF));
        }
    }

    public boolean hasData() { return !mBedRegions.isEmpty() || !mVcfEntries.isEmpty(); }
    public boolean hasValidData() { return mHasValidData; }

    public void annotateVariant(final VariantData variant)
    {
        if(mBedRegions.stream().anyMatch(x -> x.containsPosition(variant.Chromosome, variant.Position)))
        {
            variant.context().getCommonInfo().putAttribute(BLACKLIST_BED_FLAG, true);
        }

        if(mVcfEntries.stream().anyMatch(x -> x.matches(variant)))
        {
            variant.context().getCommonInfo().putAttribute(BLACKLIST_VCF_FLAG, true);
        }
    }

    private void loadBedEntries(final String filename)
    {
        if(filename == null)
            return;

        if(!Files.exists(Paths.get(filename)))
        {
            mHasValidData = false;
            return;
        }

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            for(String line : lines)
            {
                final String[] values = line.split("\t", -1);
                mBedRegions.add(new ChrBaseRegion(
                        values[0], Integer.parseInt(values[1]) + 1, Integer.parseInt(values[2])));
            }

            PV_LOGGER.info("loaded {} Blacklist BED entries from file({})", mBedRegions.size(), filename);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load Blacklist BED file: {}", e.toString());
            mHasValidData = false;
        }
    }

    private void loadVcfEntries(final String filename)
    {
        if(filename == null)
            return;

        if(!Files.exists(Paths.get(filename)))
        {
            mHasValidData = false;
            return;
        }

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                    filename, new VCFCodec(), false);

            for(VariantContext context : reader.iterator())
            {
                int position = context.getStart();
                String ref = context.getReference().getBaseString();
                String alt = context.getAlternateAlleles().get(0).toString();

                mVcfEntries.add(new VcfEntry(position, ref, alt));
            }

            PV_LOGGER.info("loaded {} BLacklist VCF entries from file({})", mVcfEntries.size(), filename);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to load Blacklist BED file: {}", e.toString());
            mHasValidData = false;
        }
    }

    public static void addHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(
                BLACKLIST_BED_FLAG, 0, VCFHeaderLineType.Flag, "Sites listed in Blacklist BED"));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                BLACKLIST_VCF_FLAG, 0, VCFHeaderLineType.Flag, "Sites listed in Blacklist VCF"));
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(BLACKLIST_BED, true, "Blacklist BED file");
        options.addOption(BLACKLIST_VCF, true, "Blacklist VCF file");
    }

    private class VcfEntry
    {
        public final int Position;
        public final String Ref;
        public final String Alt;

        public VcfEntry(final int position, final String ref, final String alt)
        {
            Position = position;
            Ref = ref;
            Alt = alt;
        }

        public boolean matches(final VariantData variant)
        {
            return variant.Position == Position && variant.Ref.equals(Ref) && variant.Alt.equals(Alt);
        }

        public String toString() { return String.format("%d %s>%s", Position, Ref, Alt); }
    }



}
