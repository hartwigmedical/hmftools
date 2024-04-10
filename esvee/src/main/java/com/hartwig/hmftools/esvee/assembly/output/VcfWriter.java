package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.NEG_ORIENTATION_CHAR;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.POS_ORIENTATION_CHAR;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.SINGLE_BREAKEND_CHAR;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ANCHOR_SUPPORT_CIGAR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ANCHOR_SUPPORT_CIGAR_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ANCHOR_SUPPORT_CIGAR_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ANCHOR_SUPPORT_CIGAR_LENGTH_DESC;
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
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REMOTE_LINKED_BY;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REMOTE_LINKED_BY_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_READS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_READS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SVTYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SVTYPE_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_FRAG_COUNT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_FRAG_COUNT_DESC;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.assembly.types.AssemblyLink;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.LinkType;
import com.hartwig.hmftools.esvee.assembly.types.PhaseSet;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;
import com.hartwig.hmftools.esvee.assembly.filters.FilterType;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class VcfWriter implements AutoCloseable
{
    private final AssemblyConfig mConfig;
    private final VariantContextWriter mWriter;

    private final List<String> mSampleNames;

    private final List<VariantContext> mVariants;

    private static final List<Allele> NO_GENOTYPE_ALLELES = Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL);

    public VcfWriter(final AssemblyConfig config)
    {
        mConfig = config;

        mSampleNames = mConfig.combinedSampleIds();
        mVariants = new ArrayList<>();

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

            writeHeader();
        }
        else
        {
            mWriter = null;
        }
    }

    private void writeHeader()
    {
        // no longer written:
        // new VCFFilterHeaderLine("GERMLINE", "This variant has support in the reference sample"),
        // new VCFInfoHeaderLine("BEID_LEN", 1, VCFHeaderLineType.Integer, "The number of assemblies associated with this variant"),

        Set<VCFHeaderLine> metaData = Sets.newHashSet();

        // TODO: check which of these are required
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
        // metaData.add(new VCFInfoHeaderLine(MAPQ, 1, VCFHeaderLineType.Integer, MAPQ_DESC));
        // metaData.add(new VCFInfoHeaderLine(BEALN, 1, VCFHeaderLineType.String, BEALN_DESC));
        //metaData.add(new VCFInfoHeaderLine(OVERHANG, 1, VCFHeaderLineType.Integer, OVERHANG_DESC));

        for(FilterType filter : FilterType.values())
        {
            metaData.add(new VCFFilterHeaderLine(filter.filterName(), filter.vcfDescription()));
        }

        // metaData.add(new VCFFormatHeaderLine(QUAL, 1, VCFHeaderLineType.Integer, QUAL_DESC));

        metaData.add(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        metaData.add(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        metaData.add(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));

        metaData.add(new VCFFormatHeaderLine(DISCORDANT_READS, 1, VCFHeaderLineType.Integer, DISCORDANT_READS_DESC));
        metaData.add(new VCFFormatHeaderLine(SPLIT_READS, 1, VCFHeaderLineType.Integer, SPLIT_READS_DESC));
        metaData.add(new VCFFormatHeaderLine(SV_FRAG_COUNT, 1, VCFHeaderLineType.Integer, SV_FRAG_COUNT_DESC));

        final VCFHeader header = new VCFHeader(metaData, mSampleNames);

        mWriter.writeHeader(header);
    }

    @Override
    public void close()
    {
        if(mWriter == null)
            return;

        // FIXME: ensure standard ordering
        // mVariants.sort(NaturalSortComparator.of(VariantContext::getContig).thenComparingInt(VariantContext::getStart));

        mVariants.forEach(mWriter::add);
        mWriter.close();
    }

    public void addVariant(final JunctionAssembly assembly)
    {
        List<SampleData> sampleDataList = buildSampleData(assembly);

        List<Genotype> genotypes = Lists.newArrayList();

        for(int i = 0; i < mSampleNames.size(); ++i)
        {
            String sampleId = mSampleNames.get(i);
            SampleData sampleData = sampleDataList.get(i);

            genotypes.add(buildGenotype(assembly, sampleId, sampleData));
        }

        PhaseSet phaseSet = assembly.phaseSet();
        List<AssemblyLink> assemblyLinks = phaseSet != null ? phaseSet.findAssemblyLinks(assembly) : Collections.emptyList();
        AssemblyLink splitLink = assemblyLinks.stream().filter(x -> x.type() == LinkType.SPLIT).findFirst().orElse(null);
        StructuralVariantType svType = splitLink != null ? splitLink.svType() : SGL;

        List<Allele> alleles = buildAlleleInfo(assembly, splitLink, svType);

        // see comments when preparing final variants and setting assembly IDs
        String assemblyId = String.valueOf(assembly.id());


        // TODO: modify post-alignment and homology sliding as required
        int breakendPosition = assembly.junction().Position;

        double mapQualityTotal = assembly.support().stream().mapToInt(x -> x.cachedRead().mappingQuality()).sum();
        int averageMapQuality = (int)Math.round(mapQualityTotal / assembly.supportCount());
        double tempQualityScore = sampleDataList.stream().mapToDouble(x -> x.SplitFragments + x.DiscordantFragments).sum();

        Set<String> filters = assembly.filters().stream().map(x -> x.filterName()).collect(Collectors.toSet());

        VariantContextBuilder builder = new VariantContextBuilder()
                .id(assemblyId)
                .chr(assembly.junction().Chromosome)
                .start(breakendPosition)
                .alleles(alleles)
                .log10PError(tempQualityScore / -10.0)
                .filters(filters)
                .genotypes(genotypes);

        mVariants.add(builder.make());

        // builder.attribute(MAP)

        /*

        builder
                .attribute(EVENT, callID)
                // .attribute("CLASSIFICATION", call.Classification.toString())
                .attribute(SVTYPE, variantCall.Classification.Type.toString())
                // .attribute("GERMLINE_SUPPORT_CNT", variantCall.germlineSupport())
                // .attribute("SOMATIC_SUPPORT_CNT", variantCall.somaticSupport())
                .attribute(MAPQ, mapQ)
                .attribute(OVERHANG, overhang)
        ;

        builder.attribute(BEID, assemblyNames);
        builder.attribute(BEIDL, left ? assemblyLeftIndices : assemblyRightIndices);
        builder.attribute(BEIDH, left ? assemblyRightIndices : assemblyLeftIndices);
        builder.attribute(ANCHOR_SUPPORT_CIGAR, left ? anchorLeftCigars : anchorRightCigars);
        builder.attribute(ANCHOR_SUPPORT_CIGAR_LENGTH, left ? anchorLeftCigarLengths : anchorRightCigarLengths);

        builder.attribute(ASSEMBLY, assemblies);
         */
    }

    private class SampleData
    {
        public int SplitFragments;
        public int DiscordantFragments;

        public SampleData()
        {
            SplitFragments = 0;
            DiscordantFragments = 0;
        }
    }

    private List<SampleData> buildSampleData(final JunctionAssembly assembly)
    {
        List<SampleData> sampleDataList = Lists.newArrayListWithExpectedSize(mSampleNames.size());

        mSampleNames.forEach(x -> sampleDataList.add(new SampleData()));

        String currentSampleId = "";
        int currentSampleIndex = 0;
        SampleData sampleData = null;

        for(SupportRead support : assembly.support())
        {
            if(support.type() == SupportType.JUNCTION_MATE)
                continue;

            int sampleIndex = support.cachedRead().sampleIndex();
            String readSampleId = mSampleNames.get(sampleIndex);

            if(!currentSampleId.equals(readSampleId))
            {
                currentSampleId = readSampleId;
                currentSampleIndex = getSampleIdIndex(readSampleId);
                sampleData = sampleDataList.get(currentSampleIndex);
            }

            if(support.type().isSplitSupport())
                ++sampleData.SplitFragments;
            else
                ++sampleData.DiscordantFragments;
        }

        return sampleDataList;
    }

    private int getSampleIdIndex(final String sampleId)
    {
        for(int i = 0; i < mSampleNames.size(); ++i)
        {
            if(mSampleNames.get(i).equals(sampleId))
                return i;
        }

        return -1;
    }

    private Genotype buildGenotype(final JunctionAssembly assembly, final String sampleId, final SampleData sampleData)
    {
        GenotypeBuilder builder = new GenotypeBuilder(sampleId);

        int depth = 0; // set from full BAM in a separate process
        int altSupport = sampleData.SplitFragments + sampleData.DiscordantFragments;

        // ideas for more per-sample data:
        // - average read quality
        // - trimmed base qualities
        // - average extension or read base distance
        // - strand bias

        builder.DP(depth)
                .AD(new int[] { 0, altSupport })
                .alleles(NO_GENOTYPE_ALLELES); // is this optional for a genotype?


        return builder.make();
    }

    private List<Allele> buildAlleleInfo(final JunctionAssembly assembly, final AssemblyLink splitLink, final StructuralVariantType svType)
    {
        // pos orientation for DEL and DUP: AGAGATTATACTTTGTGTA[10:89712341[
        // pos orientation for INV: G]3:26664499]

        // neg orientation for DEL and DUP: ]10:89700299]GAGATTATACTTTGTGTAA
        // neg orientation for INV: [3:24566181[C

        // FIXME: see Gripss method formPairedAltString() and share code if possible

        // position orientation for SG: extension sequence then '.'
        // negative orientation for SG: '.' then extension sequence
        byte[] refBase = { assembly.bases()[assembly.junctionIndex()] };
        Allele refAllele = Allele.create(refBase, true);

        String insertedBases = "";

        if(splitLink != null)
        {
            int insertLength = splitLink.insertedBases().length();

            if(insertLength > 0)
            {
                String extensionBases = assembly.formJunctionSequence();

                if(assembly.isForwardJunction())
                {
                    insertedBases = extensionBases.substring(0, insertLength);
                }
                else
                {
                    int extBaseLength = extensionBases.length();
                    insertedBases = extensionBases.substring(extBaseLength - insertLength);
                }
            }
        }
        else
        {
            insertedBases = assembly.formJunctionSequence();
        }

        StringBuilder altBases = new StringBuilder();

        if(svType != SGL)
        {
            JunctionAssembly otherAssembly = splitLink.otherAssembly(assembly);

            // TODO: see previous comment about using correct breakend position
            String otherBreakendInfo = formOtherBreakendInfo(
                    otherAssembly.junction().Chromosome, otherAssembly.junction().Position, otherAssembly.junction().Orientation);

            if(assembly.isForwardJunction())
            {
                altBases.append((char)refBase[0]);

                if(!insertedBases.isEmpty())
                    altBases.append(insertedBases);

                // details about the other breakend
                altBases.append(otherBreakendInfo);
            }
            else
            {
                altBases.append(otherBreakendInfo);

                // FIXME: check orientation or take from the extension bses
                if(!insertedBases.isEmpty())
                    altBases.append(insertedBases);

                altBases.append((char)refBase[0]);
            }
        }
        else
        {
            if(assembly.isForwardJunction())
            {
                altBases.append(insertedBases);
                altBases.append(SINGLE_BREAKEND_CHAR);
            }
            else
            {
                altBases.append(SINGLE_BREAKEND_CHAR);
                altBases.append(insertedBases);
            }
        }

        Allele altAllele = Allele.create(altBases.toString().getBytes(), false);

        return List.of(refAllele, altAllele);
    }

    private static String formOtherBreakendInfo(final String chromosome, final int position, final byte orientation)
    {
        // ]3:26664499]
        char orientationChar = orientation == POS_ORIENT ? POS_ORIENTATION_CHAR : NEG_ORIENTATION_CHAR;
        return format("%c%s:%d%c", orientationChar, chromosome, position, orientationChar);
    }
}
