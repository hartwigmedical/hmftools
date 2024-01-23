package com.hartwig.hmftools.esvee.output;

import static com.hartwig.hmftools.common.sv.SvVcfTags.ANCHOR_SUPPORT_CIGAR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ANCHOR_SUPPORT_CIGAR_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ANCHOR_SUPPORT_CIGAR_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ANCHOR_SUPPORT_CIGAR_LENGTH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASSEMBLY;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASSEMBLY_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BEALN;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BEALN_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BEID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BEIDH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BEIDH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BEIDL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BEIDL_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BEID_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISCORDANT_READS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISCORDANT_READS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.EVENT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.EVENT_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LOCAL_LINKED_BY;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LOCAL_LINKED_BY_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MAPQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MAPQ_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.OVERHANG;
import static com.hartwig.hmftools.common.sv.SvVcfTags.OVERHANG_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.QUAL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.QUAL_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REMOTE_LINKED_BY;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REMOTE_LINKED_BY_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_READS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_READS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SVTYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SVTYPE_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_FRAG_COUNT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_FRAG_COUNT_DESC;
import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.filters.SoftFilters.applyFilters;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.WriteType;
import com.hartwig.hmftools.esvee.filters.FilterType;
import com.hartwig.hmftools.esvee.old.VariantAssembly;
import com.hartwig.hmftools.esvee.util.NaturalSortComparator;
import com.hartwig.hmftools.esvee.variant.VariantCall;

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
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfWriter implements AutoCloseable
{
    private final SvConfig mConfig;
    private final VariantContextWriter mWriter;
    private final List<VariantContext> mVariants = new ArrayList<>();

    public VcfWriter(final SvConfig config)
    {
        mConfig = config;

        if(config.WriteTypes.contains(WriteType.VCF) && config.VcfFile != null)
        {
            final RefGenomeSource refGenomeSource = (RefGenomeSource) config.RefGenome;

            final SAMSequenceDictionary sequenceDictionary = refGenomeSource.refGenomeFile().getSequenceDictionary();

            mWriter = new VariantContextWriterBuilder()
                    .setOutputFile(config.VcfFile)
                    .modifyOption(Options.INDEX_ON_THE_FLY, true)
                    .modifyOption(Options.USE_ASYNC_IO, false)
                    .setReferenceDictionary(sequenceDictionary)
                    .build();

            writeHeader(mConfig.SampleNames);
        }
        else
        {
            mWriter = null;
        }
    }

    private void writeHeader(final List<String> sampleNames)
    {
        // no longer written:
        // new VCFFilterHeaderLine("GERMLINE", "This variant has support in the reference sample"),
        // new VCFFilterHeaderLine("LIKELY_FALSE", "If this variant is most likely a false-positive"),
        // new VCFInfoHeaderLine("CLASSIFICATION", 1, VCFHeaderLineType.String, ""),
        // new VCFInfoHeaderLine("GERMLINE_SUPPORT_CNT", 1, VCFHeaderLineType.Integer, "Amount of support in samples identified as germline"),
        // new VCFInfoHeaderLine("SOMATIC_SUPPORT_CNT", 1, VCFHeaderLineType.Integer, "Amount of support in samples not identified as germline"),
        // new VCFInfoHeaderLine("BEID_LEN", 1, VCFHeaderLineType.Integer, "The number of assemblies associated with this variant"),

        Set<VCFHeaderLine> metaData = Sets.newHashSet();

        metaData.add(new VCFInfoHeaderLine(EVENT, 1, VCFHeaderLineType.String, EVENT_DESC));
        metaData.add(new VCFInfoHeaderLine(MATE_ID, 1, VCFHeaderLineType.String, MATE_ID_DESC));
        metaData.add(new VCFInfoHeaderLine(SVTYPE, 1, VCFHeaderLineType.String, SVTYPE_DESC));
        metaData.add(new VCFInfoHeaderLine(BEID, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, BEID_DESC));
        metaData.add(new VCFInfoHeaderLine(BEIDL, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, BEIDL_DESC));
        metaData.add(new VCFInfoHeaderLine(BEIDH, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, BEIDH_DESC));
        metaData.add(new VCFInfoHeaderLine(ANCHOR_SUPPORT_CIGAR, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, ANCHOR_SUPPORT_CIGAR_DESC));
        metaData.add(new VCFInfoHeaderLine(ANCHOR_SUPPORT_CIGAR_LENGTH, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, ANCHOR_SUPPORT_CIGAR_LENGTH_DESC));
        metaData.add(new VCFInfoHeaderLine(LOCAL_LINKED_BY,  VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, LOCAL_LINKED_BY_DESC));
        metaData.add(new VCFInfoHeaderLine(REMOTE_LINKED_BY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, REMOTE_LINKED_BY_DESC));
        metaData.add(new VCFInfoHeaderLine(CIPOS, 2, VCFHeaderLineType.Integer, CIPOS_DESC));
        metaData.add(new VCFInfoHeaderLine(HOMSEQ, 1, VCFHeaderLineType.String, HOMSEQ_DESC));
        metaData.add(new VCFInfoHeaderLine(MAPQ, 1, VCFHeaderLineType.Integer, MAPQ_DESC));
        metaData.add(new VCFInfoHeaderLine(BEALN, 1, VCFHeaderLineType.String, BEALN_DESC));
        metaData.add(new VCFInfoHeaderLine(OVERHANG, 1, VCFHeaderLineType.Integer, OVERHANG_DESC));
        metaData.add(new VCFInfoHeaderLine(ASSEMBLY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, ASSEMBLY_DESC));

        for(FilterType filter : FilterType.values())
        {
            metaData.add(new VCFFilterHeaderLine(filter.filterName(), filter.vcfDescription()));
        }

        metaData.add(new VCFFormatHeaderLine(DISCORDANT_READS, 1, VCFHeaderLineType.Integer, DISCORDANT_READS_DESC));
        metaData.add(new VCFFormatHeaderLine(QUAL, 1, VCFHeaderLineType.Integer, QUAL_DESC));
        metaData.add(new VCFFormatHeaderLine(SPLIT_READS, 1, VCFHeaderLineType.Integer, SPLIT_READS_DESC));
        metaData.add(new VCFFormatHeaderLine(SV_FRAG_COUNT, 1, VCFHeaderLineType.Integer, SV_FRAG_COUNT_DESC));

        final VCFHeader header = new VCFHeader(metaData, sampleNames);

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
        if(mWriter == null)
            return;

        mVariants.sort(NaturalSortComparator.of(VariantContext::getContig).thenComparingInt(VariantContext::getStart));
        mVariants.forEach(mWriter::add);
        mWriter.close();
    }

    private VariantContextBuilder variant(final VariantCall variantCall, final int index, final boolean left)
    {
        final String callID = variantCall.compactName();
        final String chromosome = left ? variantCall.LeftChromosome : variantCall.RightChromosome;
        final int position = left ? variantCall.LeftPosition : variantCall.RightPosition;
        final String ref = left ? variantCall.leftRef(mConfig.RefGenome) : variantCall.rightRef(mConfig.RefGenome);
        final String descriptor = left ? variantCall.LeftDescriptor : variantCall.RightDescriptor;
        assert descriptor != null;
        final int mapQ = left ? variantCall.LeftMappingQuality : variantCall.RightMappingQuality;
        final int overhang = variantCall.overhang();

        final List<Genotype> genotypes = variantCall.sampleSupport().stream()
                .map(support -> new GenotypeBuilder(support.sampleName())
                        .DP(support.discordantPairFragmentCount())
                        .attribute(QUAL, support.quality())
                        .attribute(SPLIT_READS, support.splitReadFragmentCount())
                        .attribute(SV_FRAG_COUNT, support.totalSupportFragmentCount())
                        .make())
                .sorted(Comparator.comparing(Genotype::getSampleName))
                .collect(Collectors.toList());

        Set<String> filters = applyFilters(variantCall);

        final VariantContextBuilder builder = new VariantContextBuilder()
                .id(callID + (index != 0 ? "_" + index : ""))
                .chr(chromosome)
                .start(position)
                .alleles(ref, descriptor)
                .log10PError(-mapQ)
                .filters(filters)
                .genotypes(genotypes);

        if(index != 0)
            builder.attribute(MATE_ID, callID + "_" + (index == 1 ? 2 : 1));

        builder
                .attribute(EVENT, callID)
                // .attribute("CLASSIFICATION", call.Classification.toString())
                .attribute(SVTYPE, variantCall.Classification.Type.toString())
                // .attribute("GERMLINE_SUPPORT_CNT", variantCall.germlineSupport())
                // .attribute("SOMATIC_SUPPORT_CNT", variantCall.somaticSupport())
                .attribute(MAPQ, mapQ)
                .attribute(OVERHANG, overhang)
        ;

        final List<String> assemblyNames = new ArrayList<>();
        final List<String> assemblies = new ArrayList<>();
        final List<Integer> assemblyLeftIndices = new ArrayList<>();
        final List<Integer> assemblyRightIndices = new ArrayList<>();
        final List<String> anchorLeftCigars = new ArrayList<>();
        final List<String> anchorRightCigars = new ArrayList<>();
        final List<Integer> anchorLeftCigarLengths = new ArrayList<>();
        final List<Integer> anchorRightCigarLengths = new ArrayList<>();

        for(VariantAssembly assembly : variantCall.variantAssemblies())
        {
            assemblyNames.add(assembly.Assembly.Name);
            assemblyLeftIndices.add(assembly.LeftPosition);
            assemblyRightIndices.add(assembly.RightPosition);
            if(assembly.LeftAnchorCigar != null)
            {
                anchorLeftCigars.add(assembly.LeftAnchorCigar.toString());
                anchorLeftCigarLengths.add(assembly.LeftCigarLength);
            }
            if(assembly.RightAnchorCigar != null)
            {
                anchorRightCigars.add(assembly.RightAnchorCigar.toString());
                anchorRightCigarLengths.add(assembly.RightCigarLength);
            }
            assemblies.add(assembly.Assembly.Assembly);
        }

        builder.attribute(BEID, assemblyNames);
        builder.attribute(BEIDL, left ? assemblyLeftIndices : assemblyRightIndices);
        builder.attribute(BEIDH, left ? assemblyRightIndices : assemblyLeftIndices);
        builder.attribute(ANCHOR_SUPPORT_CIGAR, left ? anchorLeftCigars : anchorRightCigars);
        builder.attribute(ANCHOR_SUPPORT_CIGAR_LENGTH, left ? anchorLeftCigarLengths : anchorRightCigarLengths);

        builder.attribute(ASSEMBLY, assemblies);

        return builder;
    }

    /*
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
    */
}
