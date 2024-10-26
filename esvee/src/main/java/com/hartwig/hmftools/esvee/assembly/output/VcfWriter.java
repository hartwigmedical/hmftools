package com.hartwig.hmftools.esvee.assembly.output;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_N_BASE;
import static com.hartwig.hmftools.common.codon.Nucleotides.isValidDnaBase;
import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALTALN;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALTALN_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_ID_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LENGTH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_SEG_INDEX;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_SEG_INDEX_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LINKS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LINKS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.AVG_FRAG_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.AVG_FRAG_LENGTH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ASM_ORIENT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ASM_ORIENT_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ASM_POS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ASM_POS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ORIENT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.BE_ORIENT_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISC_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISC_FRAGS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ESVEE_VERSION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INSALN;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INSALN_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.MATE_ID_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_ID_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.VCF_ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.version.VersionInfo.fromAppName;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_ALIGN_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_ALIGN_LENGTH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_ID;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_ID_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_MAPQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_MAPQ_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_REPEAT_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_REPEAT_LENGTH_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_SCORE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_SCORE_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_TYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_TYPE_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS_DESC;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formPairedAltString;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formSingleAltString;
import static com.hartwig.hmftools.esvee.alignment.AlternativeAlignment.toVcfTag;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.alignment.AlternativeAlignment;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.alignment.Breakend;
import com.hartwig.hmftools.esvee.alignment.BreakendSegment;
import com.hartwig.hmftools.esvee.alignment.BreakendSupport;
import com.hartwig.hmftools.esvee.common.FilterType;

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
    private final Map<String,Integer> mSampleNameIndex; // back to config sample ordering, so genotype order can differ

    private final List<VariantContext> mVariants;

    private static final List<Allele> NO_GENOTYPE_ALLELES = Lists.newArrayList(Allele.NO_CALL, Allele.NO_CALL);

    public VcfWriter(final AssemblyConfig config)
    {
        mConfig = config;

        // reference sample IDs will be written first
        mSampleNames = Lists.newArrayList(mConfig.ReferenceIds);
        mSampleNames.addAll(mConfig.TumorIds);

        mSampleNameIndex = Maps.newHashMap();

        List<String> configSampleIds = mConfig.combinedSampleIds();
        for(int i = 0; i < configSampleIds.size(); ++i)
        {
            String sampleId = configSampleIds.get(i);
            mSampleNameIndex.put(sampleId, i);
        }

        mVariants = new ArrayList<>();

        if(config.WriteTypes.contains(WriteType.VCF))
        {
            String vcfFilename = mConfig.outputFilename(WriteType.VCF);
            final RefGenomeSource refGenomeSource = (RefGenomeSource) config.RefGenome;

            SAMSequenceDictionary sequenceDictionary = refGenomeSource.refGenomeFile().getSequenceDictionary();

            mWriter = new VariantContextWriterBuilder()
                    .setOutputFile(vcfFilename)
                    .modifyOption(Options.INDEX_ON_THE_FLY, true)
                    .modifyOption(Options.USE_ASYNC_IO, false)
                    .setReferenceDictionary(sequenceDictionary)
                    .build();

            writeHeader(sequenceDictionary);
        }
        else
        {
            mWriter = null;
        }
    }

    private void writeHeader(final SAMSequenceDictionary sequenceDictionary)
    {
        Set<VCFHeaderLine> metaData = Sets.newHashSet();

        VersionInfo versionInfo = fromAppName(APP_NAME);
        metaData.add(new VCFHeaderLine(ESVEE_VERSION, versionInfo.version()));

        metaData.add(new VCFFormatHeaderLine(QUAL, 1, VCFHeaderLineType.Integer, QUAL_DESC));

        metaData.add(new VCFInfoHeaderLine(MATE_ID, 1, VCFHeaderLineType.String, MATE_ID_DESC));
        metaData.add(new VCFInfoHeaderLine(SV_ID, 1, VCFHeaderLineType.String, SV_ID_DESC));
        metaData.add(new VCFInfoHeaderLine(SV_TYPE, 1, VCFHeaderLineType.String, SV_TYPE_DESC));

        metaData.add(new VCFInfoHeaderLine(CIPOS, 2, VCFHeaderLineType.Integer, CIPOS_DESC));
        metaData.add(new VCFInfoHeaderLine(HOMSEQ, 1, VCFHeaderLineType.String, HOMSEQ_DESC));
        metaData.add(new VCFInfoHeaderLine(IHOMPOS, 2, VCFHeaderLineType.Integer, IHOMPOS_DESC));

        metaData.add(new VCFInfoHeaderLine(INSALN, 1, VCFHeaderLineType.String, INSALN_DESC));
        metaData.add(new VCFInfoHeaderLine(ALTALN, 1, VCFHeaderLineType.String, ALTALN_DESC));
        metaData.add(new VCFInfoHeaderLine(HOMSEQ, 1, VCFHeaderLineType.String, HOMSEQ_DESC));

        metaData.add(new VCFInfoHeaderLine(SPLIT_FRAGS, 1, VCFHeaderLineType.Integer, SPLIT_FRAGS_DESC));
        metaData.add(new VCFInfoHeaderLine(DISC_FRAGS, 1, VCFHeaderLineType.Integer, DISC_FRAGS_DESC));
        metaData.add(new VCFInfoHeaderLine(TOTAL_FRAGS, 1, VCFHeaderLineType.Integer, TOTAL_FRAGS_DESC));
        metaData.add(new VCFInfoHeaderLine(AVG_FRAG_LENGTH, 1, VCFHeaderLineType.Integer, AVG_FRAG_LENGTH_DESC));
        metaData.add(new VCFInfoHeaderLine(LINE_SITE, 1, VCFHeaderLineType.Flag, LINE_SITE_DESC));

        metaData.add(new VCFInfoHeaderLine(ASM_ID, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, ASM_ID_DESC));
        metaData.add(new VCFInfoHeaderLine(ASM_LENGTH, 1, VCFHeaderLineType.Integer, ASM_LENGTH_DESC));
        metaData.add(new VCFInfoHeaderLine(ASM_SEG_INDEX, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, ASM_SEG_INDEX_DESC));
        metaData.add(new VCFInfoHeaderLine(BE_ASM_POS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, BE_ASM_POS_DESC));
        metaData.add(new VCFInfoHeaderLine(BE_ORIENT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, BE_ORIENT_DESC));
        metaData.add(new VCFInfoHeaderLine(BE_ASM_ORIENT, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, BE_ASM_ORIENT_DESC));
        metaData.add(new VCFInfoHeaderLine(SEG_ID, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, SEG_ID_DESC));
        metaData.add(new VCFInfoHeaderLine(ASM_LINKS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, ASM_LINKS_DESC));
        metaData.add(new VCFInfoHeaderLine(SEG_ALIGN_LENGTH, 1, VCFHeaderLineType.Integer, SEG_ALIGN_LENGTH_DESC));
        metaData.add(new VCFInfoHeaderLine(SEG_MAPQ, 1, VCFHeaderLineType.Integer, SEG_MAPQ_DESC));
        metaData.add(new VCFInfoHeaderLine(SEG_SCORE, 1, VCFHeaderLineType.Integer, SEG_SCORE_DESC));
        metaData.add(new VCFInfoHeaderLine(SEG_REPEAT_LENGTH, 1, VCFHeaderLineType.Integer, SEG_REPEAT_LENGTH_DESC));

        for(FilterType filter : FilterType.values())
        {
            metaData.add(new VCFFilterHeaderLine(filter.vcfTag(), filter.vcfDesc()));
        }

        metaData.add(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_KEY)));
        metaData.add(VCFStandardHeaderLines.getFormatLine((VCFConstants.GENOTYPE_ALLELE_DEPTHS)));
        metaData.add(VCFStandardHeaderLines.getFormatLine((VCFConstants.DEPTH_KEY)));

        // per sample
        metaData.add(new VCFFormatHeaderLine(SPLIT_FRAGS, 1, VCFHeaderLineType.Integer, SPLIT_FRAGS_DESC));
        metaData.add(new VCFFormatHeaderLine(DISC_FRAGS, 1, VCFHeaderLineType.Integer, DISC_FRAGS_DESC));
        metaData.add(new VCFFormatHeaderLine(TOTAL_FRAGS, 1, VCFHeaderLineType.Integer, TOTAL_FRAGS_DESC));
        metaData.add(new VCFFormatHeaderLine(STRAND_BIAS, 1, VCFHeaderLineType.Float, STRAND_BIAS_DESC));

        VCFHeader header = new VCFHeader(metaData, mSampleNames);

        header.setSequenceDictionary(sequenceDictionary);

        mWriter.writeHeader(header);
    }

    @Override
    public void close()
    {
        if(mWriter == null)
            return;

        mVariants.forEach(mWriter::add);
        mWriter.close();
    }

    private static final int[] NO_HOMOLOGY = new int[] { 0, 0 };

    public void addBreakend(final Breakend breakend)
    {
        List<Genotype> genotypes = Lists.newArrayList();

        int totalSplitFrags = 0;
        int totalDiscFrags = 0;
        for(String sampleId : mSampleNames)
        {
            int sampleSupportIndex = mSampleNameIndex.get(sampleId);
            BreakendSupport breakendSupport = breakend.sampleSupport().get(sampleSupportIndex);
            genotypes.add(buildGenotype(breakend, sampleId, breakendSupport));

            totalSplitFrags += breakendSupport.SplitFragments;
            totalDiscFrags += breakendSupport.DiscordantFragments;
        }

        int qual = breakend.calcSvQual();

        List<Allele> alleles = buildAlleleInfo(breakend);

        VariantContextBuilder builder = new VariantContextBuilder()
                .id(String.valueOf(breakend.id()))
                .chr(breakend.Chromosome)
                .start(breakend.Position)
                .alleles(alleles)
                .log10PError(qual / -10.0)
                // .filters(PASS)
                .genotypes(genotypes);

        if(!breakend.isSingle())
        {
            int otherBreakendId = breakend.otherBreakend().id();
            builder.attribute(MATE_ID, String.valueOf(otherBreakendId));

            // always write lower first
            int lowerBreakendId = otherBreakendId > breakend.id() ? breakend.id() : otherBreakendId;
            int upperBreakendId = lowerBreakendId == breakend.id() ? otherBreakendId : breakend.id();
            String svId = format("%d_%d", lowerBreakendId, upperBreakendId);

            builder.attribute(SV_ID, String.valueOf(svId));
        }

        builder.attribute(SV_TYPE, breakend.svType());

        if(breakend.Homology != null)
        {
            builder.attribute(CIPOS, new int[] { breakend.Homology.ExactStart, breakend.Homology.ExactEnd });
            builder.attribute(IHOMPOS, new int[] { breakend.Homology.InexactStart, breakend.Homology.InexactEnd });

            if(!breakend.Homology.Homology.isEmpty())
                builder.attribute(HOMSEQ, breakend.Homology.Homology);
        }
        else
        {
            // CHECK: are these optional?
            builder.attribute(CIPOS, NO_HOMOLOGY);
            builder.attribute(IHOMPOS, NO_HOMOLOGY);
        }

        builder.attribute(SPLIT_FRAGS, totalSplitFrags);
        builder.attribute(DISC_FRAGS, totalDiscFrags);
        builder.attribute(TOTAL_FRAGS, totalSplitFrags + totalDiscFrags);
        builder.attribute(AVG_FRAG_LENGTH, breakend.averageFragmentLength());

        List<AlternativeAlignment> altAlignments = breakend.alternativeAlignments();
        if(!altAlignments.isEmpty())
            builder.attribute(INSALN, toVcfTag(altAlignments));

        List<AlternativeAlignment> lowQualAltAlignments = breakend.lowQualAltAlignments();
        if(!lowQualAltAlignments.isEmpty())
            builder.attribute(ALTALN, toVcfTag(lowQualAltAlignments));

        AssemblyAlignment assemblyAlignment = breakend.assembly();

        builder.attribute(ASM_ID, assemblyAlignment.id());
        builder.attribute(ASM_LENGTH, assemblyAlignment.fullSequenceLength());

        if(assemblyAlignment.assemblies().stream().anyMatch(x -> x.hasLineSequence())
        && isMobileLineElement(breakend.Orient, breakend.InsertedBases))
        {
            builder.attribute(LINE_SITE, true);
        }

        List<BreakendSegment> segments = breakend.segments();

        builder.attribute(ASM_SEG_INDEX, segments.stream().map(x -> String.valueOf(x.Index)).collect(Collectors.joining(VCF_ITEM_DELIM)));
        builder.attribute(BE_ASM_POS, segments.stream().map(x -> String.valueOf(x.SequenceIndex)).collect(Collectors.joining(VCF_ITEM_DELIM)));
        builder.attribute(BE_ORIENT, breakend.Orient.asByte());
        builder.attribute(BE_ASM_ORIENT, segments.stream().map(x -> String.valueOf(x.Alignment.orientation().asByte())).collect(Collectors.joining(VCF_ITEM_DELIM)));

        builder.attribute(SEG_ID, segments.stream().map(x -> x.uniqueId()).collect(Collectors.joining(VCF_ITEM_DELIM)));

        // NOTE: this is used by Linx to form assembly TIs
        if(!breakend.facingBreakends().isEmpty())
        {
            builder.attribute(
                    ASM_LINKS,
                    breakend.facingBreakends().stream().map(x -> String.valueOf(x.id())).collect(Collectors.joining(VCF_ITEM_DELIM)));
        }

        builder.attribute(SEG_ALIGN_LENGTH, segments.stream().map(x -> String.valueOf(x.Alignment.alignedBases())).collect(Collectors.joining(VCF_ITEM_DELIM)));
        builder.attribute(SEG_MAPQ, segments.stream().mapToInt(x -> x.Alignment.mapQual()).max().orElse(0));
        builder.attribute(SEG_SCORE, segments.stream().mapToInt(x -> x.Alignment.score()).max().orElse(0));
        builder.attribute(SEG_REPEAT_LENGTH, segments.stream().mapToInt(x -> x.Alignment.adjustedAlignment()).max().orElse(0));

        VariantContext variantContext = builder.make();

        mVariants.add(variantContext);
    }

    private Genotype buildGenotype(final Breakend breakend, final String sampleId, final BreakendSupport breakendSupport)
    {
        GenotypeBuilder builder = new GenotypeBuilder(sampleId);

        int depth = 0; // set later by depth annotation
        int altSupport = breakendSupport.totalSupport();

        builder.DP(depth)
                .AD(new int[] { 0, altSupport }) // set later by depth annotation
                .alleles(NO_GENOTYPE_ALLELES);

        builder.attribute(SPLIT_FRAGS, breakendSupport.SplitFragments);
        builder.attribute(DISC_FRAGS, breakendSupport.DiscordantFragments);
        builder.attribute(TOTAL_FRAGS, altSupport);
        builder.attribute(STRAND_BIAS, breakendSupport.strandBias());

        return builder.make();
    }

    private List<Allele> buildAlleleInfo(final Breakend breakend)
    {
        byte[] refBases = mConfig.RefGenome.getBases(breakend.Chromosome, breakend.Position, breakend.Position);
        byte refBase = refBases[0];

        if(!isValidDnaBase(refBase))
            refBase = DNA_N_BASE;

        Allele refAllele = Allele.create(refBase, true);

        String altBase = String.valueOf((char)refBase);

        String altBases;

        if(breakend.isSingle())
        {
            altBases = formSingleAltString(altBase, breakend.InsertedBases, breakend.Orient);
        }
        else
        {
            final Breakend otherBreakend = breakend.otherBreakend();

            altBases = formPairedAltString(
                    altBase, breakend.InsertedBases, otherBreakend.Chromosome, otherBreakend.Position, breakend.Orient, otherBreakend.Orient);
        }

        Allele altAllele = Allele.create(altBases.toString().getBytes(), false);

        return List.of(refAllele, altAllele);
    }
}
