package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;
import static htsjdk.variant.vcf.VCFHeaderLineCount.UNBOUNDED;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;

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

public class ClinvarAnnotation
{
    private final Map<String,List<ClinvarEntry>> mChrEntries;
    private boolean mHasValidData;

    private static final String CLINVAR_VCF = "clinvar_vcf";

    public static final String CLNSIG = "CLNSIG";
    public static final String CLNSIGCONF = "CLNSIGCONF";

    public static final String CLNSIG_DESC = "Clinical significance for this single variant";
    public static final String CLNSIGCONF_DESC = "Conflicting clinical significance for this single variant";

    public ClinvarAnnotation(final CommandLine cmd)
    {
        mChrEntries = Maps.newHashMap();
        mHasValidData = true;

        if(cmd.hasOption(CLINVAR_VCF))
        {
            loadEntries(cmd.getOptionValue(CLINVAR_VCF));
        }
    }

    public boolean hasData() { return !mChrEntries.isEmpty(); }
    public boolean hasValidData() { return mHasValidData; }

    public void annotateVariant(final VariantData variant)
    {
        String chromosome = RefGenomeFunctions.stripChrPrefix(variant.Chromosome);

        List<ClinvarEntry> entries = mChrEntries.get(chromosome);

        if(entries == null)
            return;

        for(ClinvarEntry entry : entries)
        {
            if(entry.matches(variant))
            {
                variant.context().getCommonInfo().putAttribute(CLNSIG, entry.Significance);

                if(!entry.Conflict.isEmpty())
                    variant.context().getCommonInfo().putAttribute(CLNSIGCONF, entry.Conflict);

                return;
            }
        }
    }

    public static void addHeader(final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(CLNSIG, UNBOUNDED, VCFHeaderLineType.String, CLNSIG_DESC));
        header.addMetaDataLine(new VCFInfoHeaderLine(CLNSIGCONF, UNBOUNDED, VCFHeaderLineType.String, CLNSIGCONF_DESC));
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(CLINVAR_VCF, true, "Clinvar annotation VCF");
    }

    private void loadEntries(final String filename)
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
                if(context.getAlleles().size() < 2)
                    continue;

                String chromosome = RefGenomeFunctions.stripChrPrefix(context.getContig());

                List<ClinvarEntry> entries = mChrEntries.get(chromosome);

                if(entries == null)
                {
                    entries = Lists.newArrayList();
                    mChrEntries.put(chromosome, entries);
                }

                int position = context.getStart();
                String ref = context.getReference().getBaseString();
                String alt = context.getAlternateAlleles().get(0).toString();

                String significance = context.getAttributeAsString(CLNSIG, "");
                String conflict = context.getAttributeAsString(CLNSIGCONF, "");

                if(significance.isEmpty() && conflict.isEmpty())
                    continue;

                entries.add(new ClinvarEntry(position, ref, alt, stripBrackets(significance), stripBrackets(conflict)));
            }

            PV_LOGGER.info("loaded {} Clinvar entries from file({})",
                    mChrEntries.values().stream().mapToInt(x -> x.size()).sum(), filename);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to read Clinvar VCF file: {}",  e.toString());
            mHasValidData = false;
        }
    }

    private static String stripBrackets(final String clinvarStr)
    {
        return clinvarStr.replaceAll("\\[", "").replaceAll("\\]", "").replaceAll(" ", "");
    }

    private class ClinvarEntry
    {
        public final int Position;
        public final String Ref;
        public final String Alt;
        public final String Significance;
        public final String Conflict;

        public ClinvarEntry(final int position, final String ref, final String alt, final String significance, final String conflict)
        {
            Position = position;
            Ref = ref;
            Alt = alt;
            Significance = significance;
            Conflict = conflict;
        }

        public boolean matches(final VariantData variant)
        {
            return variant.Position == Position && variant.Ref.equals(Ref) && variant.Alt.equals(Alt);
        }

        public String toString() { return String.format("%d %s>%s details(%s - %s)",
                Position, Ref, Alt, Significance, Conflict); }
    }
}
