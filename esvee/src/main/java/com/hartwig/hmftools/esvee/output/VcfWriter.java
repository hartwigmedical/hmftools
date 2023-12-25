package com.hartwig.hmftools.esvee.output;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.Context;
import com.hartwig.hmftools.esvee.SvConstants;
import com.hartwig.hmftools.esvee.models.AssemblyClassification;
import com.hartwig.hmftools.esvee.util.NaturalSortComparator;
import com.hartwig.hmftools.esvee.processor.VariantCall;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfWriter implements AutoCloseable
{
    private final Context mContext;
    private final VariantContextWriter mWriter;
    private final List<VariantContext> mVariants = new ArrayList<>();

    public VcfWriter(final Context context, final List<String> sampleNames)
    {
        mContext = context;

        final SAMSequenceDictionary sequenceDictionary = context.ReferenceGenome.refGenomeFile().getSequenceDictionary();
        mWriter = new VariantContextWriterBuilder()
                .setOutputFile(context.Config.VcfFile)
                .modifyOption(Options.INDEX_ON_THE_FLY, false)
                .modifyOption(Options.USE_ASYNC_IO, false)
                .setReferenceDictionary(sequenceDictionary)
                .build();

        final VCFHeader header = new VCFHeader(new LinkedHashSet<>(List.of(
                new VCFFilterHeaderLine("GERMLINE", "This variant has support in the reference sample"),
                new VCFFilterHeaderLine("MULTIPLE_ASSEMBLIES", "Whether this variant is supported by multiple assemblies"),
                new VCFFilterHeaderLine("LOW_OVERHANG", "If split-reads do not cover this variant sufficiently"),
                new VCFFilterHeaderLine("LOW_QUALITY", "If this variant has a poor quality metric"),
                new VCFFilterHeaderLine("LOW_SUPPORT", "If this variant has low support"),
                new VCFFilterHeaderLine("LIKELY_FALSE", "If this variant is most likely a false-positive"),
                new VCFInfoHeaderLine("EVENT", 1, VCFHeaderLineType.String, ""),
                new VCFInfoHeaderLine("MATEID", 1, VCFHeaderLineType.String, ""),
                new VCFInfoHeaderLine("CLASSIFICATION", 1, VCFHeaderLineType.String, ""),
                new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String, "Variant Type"),
                new VCFInfoHeaderLine("GERMLINE_SUPPORT_CNT", 1, VCFHeaderLineType.Integer, "Amount of support in samples identified as germline"),
                new VCFInfoHeaderLine("SOMATIC_SUPPORT_CNT", 1, VCFHeaderLineType.Integer, "Amount of support in samples not identified as germline"),
                new VCFInfoHeaderLine("BEID", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Assembly Names"),
                new VCFInfoHeaderLine("BEID_LEN", 1, VCFHeaderLineType.Integer, "The number of assemblies associated with this variant"),
                new VCFInfoHeaderLine("BEIDL", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Local offset in assemblies (one per assembly)"),
                new VCFInfoHeaderLine("BEIDH", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Remote offset in assemblies (one per assembly)"),
                new VCFInfoHeaderLine("SC", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "CIGAR of local anchor (one per assembly)"),
                new VCFInfoHeaderLine("SC_LEN", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "Length of local anchor CIGAR (one per assembly)"),
                new VCFInfoHeaderLine("LOCAL_LINKED_BY",  VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "ID of previous phased variant (one per assembly)"),
                new VCFInfoHeaderLine("REMOTE_LINKED_BY", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "ID of next phased variant (one per assembly)"),
                new VCFInfoHeaderLine("CIPOS", 2, VCFHeaderLineType.Integer, "If homology is present, the uncertainty window"),
                new VCFInfoHeaderLine("HOMSEQ", 1, VCFHeaderLineType.String, "If homology is present, the bases contained in the uncertainty window"),
                new VCFInfoHeaderLine("MAPQ", 1, VCFHeaderLineType.Integer, "The anchor's mapping quality"),
                new VCFInfoHeaderLine("BEALN", 1, VCFHeaderLineType.String, "For singles, the potential alignments of the insert sequence"),
                new VCFInfoHeaderLine("OVERHANG", 1, VCFHeaderLineType.Integer, "The minimum of the left & right overhang for this assembly"),
                new VCFInfoHeaderLine("ASSEMBLY", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "The minimum of the left & right overhang for this assembly"),
                new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "The number of discordant fragments supporting this variant"),
                new VCFFormatHeaderLine("QUAL", 1, VCFHeaderLineType.Integer, ""),
                new VCFFormatHeaderLine("SR", 1, VCFHeaderLineType.Integer, "The number of split-read fragments supporting this variant"),
                new VCFFormatHeaderLine("VF", 1, VCFHeaderLineType.Integer, "The sum of SR and DP")
        )), new LinkedHashSet<>(sampleNames));

        mWriter.writeHeader(header);
    }

    public void append(final VariantCall call)
    {
        if(!call.isSingleSided())
        {
            final VariantContextBuilder left = variant(call, 1, true);
            final VariantContextBuilder right = variant(call, 2, false);
            mVariants.add(left.make());
            mVariants.add(right.make());
        }
        else
        {
            if(call.LeftDescriptor == null)
            {
                SV_LOGGER.error("Expected descriptor for {}", call);
                return;
            }
            final VariantContextBuilder left = variant(call, 0, true);
            mVariants.add(left.make());
        }
    }

    @Override
    public void close()
    {
        mVariants.sort(NaturalSortComparator.of(VariantContext::getContig).thenComparingInt(VariantContext::getStart));
        mVariants.forEach(mWriter::add);
        mWriter.close();
    }

    private VariantContextBuilder variant(final VariantCall call, final int index, final boolean left)
    {
        final String callID = call.compactName();
        final String chromosome = left ? call.LeftChromosome : call.RightChromosome;
        final int position = left ? call.LeftPosition : call.RightPosition;
        final String ref = left ? call.leftRef(mContext) : call.rightRef(mContext);
        final String descriptor = left ? call.LeftDescriptor : call.RightDescriptor;
        assert descriptor != null;
        final int mapQ = left ? call.LeftMappingQuality : call.RightMappingQuality;
        final int overhang = call.overhang();

        final List<Genotype> genotypes = call.sampleSupport().stream()
                .map(support -> new GenotypeBuilder(support.sampleName())
                        .DP(support.discordantPairFragmentCount())
                        .attribute("QUAL", support.quality())
                        .attribute("SR", support.splitReadFragmentCount())
                        .attribute("VF", support.totalSupportFragmentCount())
                        .make())
                .sorted(Comparator.comparing(Genotype::getSampleName))
                .collect(Collectors.toList());

        final Set<String> filters = new HashSet<>();
        final boolean isLowOverhang = overhang < SvConstants.VCFLOWOVERHANGTHRESHOLD;
        final boolean isLowQuality = call.quality() < SvConstants.VCFLOWQUALITYTHRESHOLD;
        final boolean isLowSupport = call.supportingFragments().size() < SvConstants.MINREADSTOSUPPORTASSEMBLY;
        final boolean isLikelyFalse = isLowSupport || (isLowOverhang && call.discordantSupport() == 0) || isLowQuality;
        if (call.isGermline())
            filters.add("GERMLINE");
        if (call.associatedAssemblies().size() > 1)
            filters.add("MULTIPLE_ASSEMBLIES");
        if (isLowOverhang)
            filters.add("LOW_OVERHANG");
        if (isLowQuality)
            filters.add("LOW_QUALITY");
        if (isLowSupport)
            filters.add("LOW_SUPPORT");
        if (isLikelyFalse)
            filters.add("LIKELY_FALSE");

        final VariantContextBuilder builder = new VariantContextBuilder()
                .id(callID + (index != 0 ? "_" + index : ""))
                .chr(chromosome)
                .start(position)
                .alleles(ref, descriptor)
                .log10PError(-mapQ)
                .filters(filters)
                .genotypes(genotypes);

        if(index != 0)
            builder.attribute("MATEID", callID + "_" + (index == 1 ? 2 : 1));

        builder
                .attribute("EVENT", callID)
                .attribute("CLASSIFICATION", call.Classification.toString())
                .attribute("SVTYPE", svType(call.Classification))
                .attribute("GERMLINE_SUPPORT_CNT", call.germlineSupport())
                .attribute("SOMATIC_SUPPORT_CNT", call.somaticSupport())
                .attribute("MAPQ", mapQ)
                .attribute("OVERHANG", overhang)
        ;

        final List<String> assemblyNames = new ArrayList<>();
        final List<String> assemblies = new ArrayList<>();
        final List<Integer> assemblyLeftIndices = new ArrayList<>();
        final List<Integer> assemblyRightIndices = new ArrayList<>();
        final List<String> anchorLeftCigars = new ArrayList<>();
        final List<String> anchorRightCigars = new ArrayList<>();
        final List<Integer> anchorLeftCigarLengths = new ArrayList<>();
        final List<Integer> anchorRightCigarLengths = new ArrayList<>();

        for (final VariantCall.VariantAssembly assembly : call.variantAssemblies())
        {
            assemblyNames.add(assembly.Assembly.Name);
            assemblyLeftIndices.add(assembly.LeftPosition);
            assemblyRightIndices.add(assembly.RightPosition);
            if (assembly.LeftAnchorCigar != null)
            {
                anchorLeftCigars.add(assembly.LeftAnchorCigar.toString());
                anchorLeftCigarLengths.add(assembly.LeftCigarLength);
            }
            if (assembly.RightAnchorCigar != null)
            {
                anchorRightCigars.add(assembly.RightAnchorCigar.toString());
                anchorRightCigarLengths.add(assembly.RightCigarLength);
            }
            assemblies.add(assembly.Assembly.Assembly);
        }

        builder.attribute("BEID", assemblyNames);
        builder.attribute("BEID_LEN", assemblyNames.size());
        builder.attribute("BEIDL", left ? assemblyLeftIndices : assemblyRightIndices);
        builder.attribute("BEIDH", left ? assemblyRightIndices : assemblyLeftIndices);
        builder.attribute("SC", left ? anchorLeftCigars : anchorRightCigars);
        builder.attribute("SC_LEN", left ? anchorLeftCigarLengths : anchorRightCigarLengths);

        // FIXME-CHASHA:
        builder.attribute("ASSEMBLY", assemblies);
        // builder.attribute("ASSEMBLY", assemblies);

        return builder;
    }

    private String svType(final AssemblyClassification classification)
    {
        switch(classification.Type)
        {
            case INSERT:
                return "INS";
            case DELETION:
                return "DEL";
            case DUPLICATION:
                return "DUP";
            case TRANSLOCATION:
                return "BND";
            case INVERSION:
                return "INV";
            case UNKNOWN:
                return "SGL";
        }
        throw new IllegalStateException("Invalid AssemblyClassification: " + classification);
    }
}
