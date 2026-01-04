package com.hartwig.hmftools.purple.sv;

import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_COUNT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_FILTER_PON;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SV_TYPE;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.TINC_RECOVERED_DESC;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.TINC_RECOVERED_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF_DESC;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_CHANGE;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_CHANGE_DESC;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_DESC;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_JUNCTION_COPY_NUMBER;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_JUNCTION_COPY_NUMBER_DESC;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import static htsjdk.variant.vcf.VCFHeaderLineCount.UNBOUNDED;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.Multimaps;
import com.hartwig.hmftools.purple.copynumber.CopyNumberEnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.apache.logging.log4j.util.Strings;

import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

public class SomaticSvCache
{
    private static final Allele REF_ALLELE = Allele.create("N", true);
    private static final Allele INCREASING_ALLELE = Allele.create(".N", false);
    private static final Allele DECREASING_ALLELE = Allele.create("N.", false);

    private final String mOutputVcfFilename;
    private final VCFHeader mVcfHeader;
    private final GenotypeIds mGenotypeIds;
    private final VariantContextCollection mVariantCollection;

    private int mNextVarId = 0;

    public SomaticSvCache()
    {
        mVcfHeader = null;
        mOutputVcfFilename = Strings.EMPTY;
        mVariantCollection = new VariantContextCollection(null);
        mGenotypeIds = null;
    }

    public SomaticSvCache(final String version, final String inputVcf, final String outputVcf, final PurpleConfig config)
    {
        VCFFileReader vcfReader = new VCFFileReader(new File(inputVcf), false);
        mOutputVcfFilename = outputVcf;
        mVcfHeader = generateOutputHeader(version, vcfReader.getFileHeader());
        mVariantCollection = new VariantContextCollection(mVcfHeader);

        mGenotypeIds = GenotypeIds.fromVcfHeader(mVcfHeader, config.ReferenceId, config.TumorId);

        for(VariantContext context : vcfReader)
        {
            mVariantCollection.addVariant(context);
        }

        vcfReader.close();

        PPL_LOGGER.info("loaded {} somatic SVs from {}", somaticVariants().size(), inputVcf);
    }

    public void inferMissingVariant(final List<PurpleCopyNumber> copyNumbers)
    {
        if(enabled())
        {
            for(int i = 1; i < copyNumbers.size(); i++)
            {
                PurpleCopyNumber copyNumber = copyNumbers.get(i);
                if(copyNumber.segmentStartSupport() == SegmentSupport.NONE)
                {
                    final PurpleCopyNumber prev = copyNumbers.get(i - 1);
                    mVariantCollection.addVariant(infer(copyNumber, prev));
                }
            }
        }
    }

    public void addTincVariant(final StructuralVariant variant)
    {
        mVariantCollection.addVariant(variant.startContext());
        variant.startContext().getCommonInfo().putAttribute(TINC_RECOVERED_FLAG, true);

        // check if should be PON filtered
        int ponCount = variant.startContext().getAttributeAsInt(PON_COUNT, 0);

        if(ponCount > 0)
            variant.startContext().getCommonInfo().addFilter(PON_FILTER_PON);

        if(variant.endContext() != null)
        {
            mVariantCollection.addVariant(variant.endContext());
            variant.endContext().getCommonInfo().putAttribute(TINC_RECOVERED_FLAG, true);

            if(ponCount > 0)
                variant.endContext().getCommonInfo().addFilter(PON_FILTER_PON);
        }
    }

    private VariantContext infer(final PurpleCopyNumber copyNumber, final PurpleCopyNumber prev)
    {
        final long position;
        final Allele allele;
        if(Doubles.greaterThan(copyNumber.averageTumorCopyNumber(), prev.averageTumorCopyNumber()))
        {
            allele = INCREASING_ALLELE;
            position = copyNumber.start();
        }
        else
        {
            allele = DECREASING_ALLELE;
            position = copyNumber.start() - 1;
        }

        final Collection<Allele> alleles = Lists.newArrayList(REF_ALLELE, allele);

        long lowerRange = Math.min(-500, copyNumber.minStart() - copyNumber.start());
        long upperRange = Math.max(500, copyNumber.maxStart() - copyNumber.start());

        return new VariantContextBuilder("purple", copyNumber.chromosome(), position, copyNumber.start(), alleles).filter(INFERRED)
                .id("purple_" + mNextVarId++)
                .attribute(CIPOS, Lists.newArrayList(lowerRange, upperRange))
                .attribute(SV_TYPE, StructuralVariantType.BND.toString())
                .attribute(INFERRED, true)
                .noGenotypes()
                .make();
    }

    public void write(
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, boolean passOnly, final Gender gender)
    {
        if(mVcfHeader == null)
            return;

        try
        {
            final VariantContextWriter writer = new VariantContextWriterBuilder()
                    .setOutputFile(mOutputVcfFilename)
                    .setReferenceDictionary(mVcfHeader.getSequenceDictionary())
                    .setIndexCreator(new TabixIndexCreator(mVcfHeader.getSequenceDictionary(), new TabixFormat()))
                    .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                    .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();

            writer.writeHeader(mVcfHeader);

            // builds a new collection with annotated variant contexts - leaves the original one as-is
            VariantContextCollection enrichedCollection = getEnrichedCollection(purityAdjuster, copyNumbers, gender);

            Iterator<VariantContext> variantIter = enrichedCollection.iterator();

            while(variantIter.hasNext())
            {
                VariantContext variant = variantIter.next();

                if(passOnly && variant.isFiltered())
                    continue;

                writer.add(variant);
            }

            writer.close();
        }
        catch(Exception e)
        {
              PPL_LOGGER.error("failed to write SV VCF({}): {}", mOutputVcfFilename, e.toString());
        }
    }

    private VariantContextCollection getEnrichedCollection(
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, final Gender gender)
    {
        StructuralVariantFactory svFactory = StructuralVariantFactory.build();
        svFactory.setGenotypeOrdinals(mGenotypeIds.ReferenceOrdinal, mGenotypeIds.TumorOrdinal);

        Iterator<VariantContext> variantIter = mVariantCollection.iterator();

        while(variantIter.hasNext())
        {
            VariantContext variant = variantIter.next();
            svFactory.addVariantContext(variant);
        }

        CopyNumberEnrichedStructuralVariantFactory svEnricher = new CopyNumberEnrichedStructuralVariantFactory(
                purityAdjuster, Multimaps.fromRegions(copyNumbers));

        VariantContextCollection enrichedCollection = new VariantContextCollection(mVcfHeader);

        List<EnrichedStructuralVariant> enrichedVariants = svEnricher.enrich(svFactory.results());

        for(EnrichedStructuralVariant enrichedSV : enrichedVariants)
        {
            if(gender == Gender.FEMALE)
            {
                if(HumanChromosome._Y.matches(enrichedSV.chromosome(true))
                || (enrichedSV.end() != null && HumanChromosome._Y.matches(enrichedSV.chromosome(false))))
                {
                    continue;
                }
            }

            addEnrichedVariantContexts(enrichedCollection, enrichedSV);
        }

        svFactory.unmatched().forEach(enrichedCollection::addVariant);

        return enrichedCollection;
    }

    public static void addEnrichedVariantContexts(final VariantContextCollection variantCollection, final EnrichedStructuralVariant variant)
    {
        final VariantContext startContext = variant.startContext();
        final VariantContext endContext = variant.endContext();

        List<Double> purpleAF = Lists.newArrayList();
        List<Double> purpleCN = Lists.newArrayList();
        List<Double> purpleCNChange = Lists.newArrayList();

        if(variant.start().adjustedAlleleFrequency() != null)
            purpleAF.add(variant.start().adjustedAlleleFrequency());

        if(variant.start().adjustedCopyNumber() != null)
            purpleCN.add(variant.start().adjustedCopyNumber());

        if(variant.start().adjustedCopyNumberChange() != null)
            purpleCNChange.add(variant.start().adjustedCopyNumberChange());

        if(variant.end() != null)
        {
            if(variant.end().adjustedAlleleFrequency() != null)
                purpleAF.add(variant.end().adjustedAlleleFrequency());

            if(variant.end().adjustedCopyNumber() != null)
                purpleCN.add(variant.end().adjustedCopyNumber());

            if(variant.end().adjustedCopyNumberChange() != null)
                purpleCNChange.add(variant.end().adjustedCopyNumberChange());
        }

        if(startContext != null)
        {
            variantCollection.addVariant(buildVariantContext(startContext, purpleAF, purpleCN, purpleCNChange, variant.junctionCopyNumber()));
        }

        if(endContext != null)
        {
            purpleAF = Lists.newArrayList(purpleAF);
            purpleCN = Lists.newArrayList(purpleCN);
            purpleCNChange = Lists.newArrayList(purpleCNChange);

            Collections.reverse(purpleAF);
            Collections.reverse(purpleCN);
            Collections.reverse(purpleCNChange);

            variantCollection.addVariant(buildVariantContext(endContext, purpleAF, purpleCN, purpleCNChange, variant.junctionCopyNumber()));
        }
    }

    public static VariantContext buildVariantContext(
            final VariantContext variantContext, final List<Double> purpleAF, final List<Double> purpleCN,
            final List<Double> purpleCNChange, final Double junctionCopyNumber)
    {
        VariantContextBuilder builder = new VariantContextBuilder(variantContext);

        if(!purpleAF.isEmpty())
            builder.attribute(PURPLE_AF, purpleAF);

        if(!purpleCN.isEmpty())
            builder.attribute(PURPLE_CN, purpleCN);

        if(!purpleCNChange.isEmpty())
            builder.attribute(PURPLE_CN_CHANGE, purpleCNChange);

        if(junctionCopyNumber != null)
            builder.attribute(PURPLE_JUNCTION_COPY_NUMBER, junctionCopyNumber);

        return builder.make();
    }

    public List<StructuralVariant> somaticVariants() { return mVariantCollection.variants(); }

    public int passingBnd()
    {
        return (int) mVariantCollection.variants()
                .stream()
                .filter(x -> x.filter() == null || Objects.equals(x.filter(), PASS_FILTER))
                .filter(x -> !x.type().equals(StructuralVariantType.INF))
                .filter(x -> !x.type().equals(StructuralVariantType.SGL))
                .count();
    }

    public static VCFHeader generateOutputHeader(final String purpleVersion, final VCFHeader template)
    {
        final VCFHeader outputVCFHeader = new VCFHeader(template.getMetaDataInInputOrder(), template.getGenotypeSamples());
        outputVCFHeader.addMetaDataLine(new VCFHeaderLine("purpleVersion", purpleVersion));

        outputVCFHeader.addMetaDataLine(VCFStandardHeaderLines.getFormatLine("GT"));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(INFERRED, 0, VCFHeaderLineType.Flag, INFERRED_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(TINC_RECOVERED_FLAG, 0, VCFHeaderLineType.Flag, TINC_RECOVERED_DESC));
        outputVCFHeader.addMetaDataLine(new VCFFilterHeaderLine(INFERRED, INFERRED_DESC));

        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_AF, UNBOUNDED, VCFHeaderLineType.Float, PURPLE_AF_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN, UNBOUNDED, VCFHeaderLineType.Float, PURPLE_CN_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_JUNCTION_COPY_NUMBER, 1, VCFHeaderLineType.Float,
                PURPLE_JUNCTION_COPY_NUMBER_DESC));
        outputVCFHeader.addMetaDataLine(new VCFInfoHeaderLine(PURPLE_CN_CHANGE,
                UNBOUNDED,
                VCFHeaderLineType.Float,
                PURPLE_CN_CHANGE_DESC));

        return outputVCFHeader;
    }

    private boolean enabled() { return mVcfHeader != null; }
}
