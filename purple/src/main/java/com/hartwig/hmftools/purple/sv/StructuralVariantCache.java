package com.hartwig.hmftools.purple.sv;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.CIPOS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.SVTYPE;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF_INFO;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_INFO;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.Optional;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.purple.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.purple.copynumber.CopyNumberEnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantHeader;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.purple.config.ReferenceData;

import org.apache.logging.log4j.util.Strings;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class StructuralVariantCache
{
    private static final Allele REF_ALLELE = Allele.create("N", true);
    private static final Allele INCREASING_ALLELE = Allele.create(".N", false);
    private static final Allele DECREASING_ALLELE = Allele.create("N.", false);

    private final String mOutputVcfFilename;
    private final Optional<VCFHeader> mVcfHeader;
    private final VariantContextCollection mVariantCollection;
    private final IndexedFastaSequenceFile mRefGenomeFile;

    private int mNextVarId = 0;

    public StructuralVariantCache()
    {
        mVcfHeader = Optional.empty();
        mOutputVcfFilename = Strings.EMPTY;
        mVariantCollection = new VariantContextCollection(null);
        mRefGenomeFile = null;
    }

    public StructuralVariantCache(
            final String version, final String templateVCF, final String outputVCF, final ReferenceData referenceData)
    {
        final VCFFileReader vcfReader = new VCFFileReader(new File(templateVCF), false);
        mOutputVcfFilename = outputVCF;
        mVcfHeader = Optional.of(generateOutputHeader(version, vcfReader.getFileHeader()));
        mVariantCollection = new VariantContextCollection(mVcfHeader.get());
        mRefGenomeFile = referenceData.RefGenome;

        for(VariantContext context : vcfReader)
        {
            mVariantCollection.add(context);
        }

        vcfReader.close();
    }

    public static StructuralVariantCache createStructuralVariantCache(
            final String inputVcf, final String outputVcf, final String purpleVersion, final ReferenceData referenceData)
    {
        if(inputVcf.isEmpty())
            return new StructuralVariantCache();

        StructuralVariantCache svCache = new StructuralVariantCache(purpleVersion, inputVcf, outputVcf, referenceData);

        PPL_LOGGER.info("loaded {} structural variants from {}", svCache.variants().size(), inputVcf);

        return svCache;
    }

    public void addVariant(final VariantContext variantContext)
    {
        if(enabled())
        {
            mVariantCollection.add(variantContext);
        }
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
                    mVariantCollection.add(infer(copyNumber, prev));
                }
            }
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
                .attribute(StructuralVariantFactory.IMPRECISE, true)
                .id("purple_" + mNextVarId++)
                .attribute(CIPOS, Lists.newArrayList(lowerRange, upperRange))
                .attribute(SVTYPE, "BND")
                .attribute(INFERRED, true)
                .noGenotypes()
                .make();
    }

    public void write(final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, boolean passOnly)
    {
        if(!mVcfHeader.isPresent())
            return;

        try
        {
            final VariantContextWriter writer = new VariantContextWriterBuilder()
                    .setOutputFile(mOutputVcfFilename)
                    .setReferenceDictionary(mVcfHeader.get().getSequenceDictionary())
                    .setIndexCreator(new TabixIndexCreator(mVcfHeader.get().getSequenceDictionary(), new TabixFormat()))
                    .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                    .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();

            final StructuralRefContextEnrichment refEnricher = new StructuralRefContextEnrichment(mRefGenomeFile, writer::add);

            writer.writeHeader(refEnricher.enrichHeader(mVcfHeader.get()));

            VariantContextCollection enrichedCollection = getEnrichedCollection(purityAdjuster, copyNumbers);

            Iterator<VariantContext> variantIter = enrichedCollection.iterator();

            while(variantIter.hasNext())
            {
                VariantContext variant = variantIter.next();

                if(passOnly && variant.isFiltered())
                    continue;

                refEnricher.accept(variant);
            }

            refEnricher.flush();

            writer.close();
        }
        catch(Exception e)
        {
              PPL_LOGGER.error("failed to write SV VCF({}): {}", mOutputVcfFilename, e.toString());
        }
    }

    private VariantContextCollection getEnrichedCollection(final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers)
    {
        final StructuralVariantFactory svFactory = new StructuralVariantFactory(x -> true);

        Iterator<VariantContext> variantIter = mVariantCollection.iterator();

        while(variantIter.hasNext())
        {
            VariantContext variant = variantIter.next();
            svFactory.addVariantContext(variant);
        }

        final CopyNumberEnrichedStructuralVariantFactory svEnricher = new CopyNumberEnrichedStructuralVariantFactory(
                purityAdjuster, Multimaps.fromRegions(copyNumbers));

        VariantContextCollection enrichedCollection = new VariantContextCollection(mVcfHeader.get());

        List<EnrichedStructuralVariant> enrichedVariants = svEnricher.enrich(svFactory.results());

        for(EnrichedStructuralVariant enrichedSV : enrichedVariants)
        {
            addEnrichedVariantContexts(enrichedCollection, enrichedSV);
        }

        svFactory.unmatched().forEach(enrichedCollection::add);

        return enrichedCollection;
    }

    private void addEnrichedVariantContexts(final VariantContextCollection enrichedCollection, final EnrichedStructuralVariant variant)
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
            enrichedCollection.add(buildVariantContext(startContext, purpleAF, purpleCN, purpleCNChange, variant.junctionCopyNumber()));
        }

        if(endContext != null)
        {
            purpleAF = Lists.newArrayList(purpleAF);
            purpleCN = Lists.newArrayList(purpleCN);
            purpleCNChange = Lists.newArrayList(purpleCNChange);

            Collections.reverse(purpleAF);
            Collections.reverse(purpleCN);
            Collections.reverse(purpleCNChange);

            enrichedCollection.add(buildVariantContext(endContext, purpleAF, purpleCN, purpleCNChange, variant.junctionCopyNumber()));
        }
    }

    private VariantContext buildVariantContext(
            final VariantContext variantContext, final List<Double> purpleAF, final List<Double> purpleCN,
            final List<Double> purpleCNChange, final Double junctionCopyNumber)
    {
        VariantContextBuilder builder = new VariantContextBuilder(variantContext);

        if(!purpleAF.isEmpty())
            builder.attribute(PURPLE_AF_INFO, purpleAF);

        if(!purpleCN.isEmpty())
            builder.attribute(PURPLE_CN_INFO, purpleCN);

        if(!purpleCNChange.isEmpty())
            builder.attribute(StructuralVariantHeader.PURPLE_CN_CHANGE_INFO, purpleCNChange);

        if(junctionCopyNumber != null)
            builder.attribute(StructuralVariantHeader.PURPLE_JUNCTION_COPY_NUMBER_INFO, junctionCopyNumber);

        return builder.make();
    }

    public List<StructuralVariant> variants() { return mVariantCollection.segmentationVariants(); }

    public int passingBnd()
    {
        return (int) mVariantCollection.segmentationVariants()
                .stream()
                .filter(x -> x.filter() == null || Objects.equals(x.filter(), PASS))
                .filter(x -> !x.type().equals(StructuralVariantType.INF))
                .filter(x -> !x.type().equals(StructuralVariantType.SGL))
                .count();
    }

    @VisibleForTesting
    static VCFHeader generateOutputHeader(final String purpleVersion, final VCFHeader template)
    {
        return StructuralVariantHeader.generateHeader(purpleVersion, template);
    }

    private boolean enabled()
    {
        return mVcfHeader.isPresent();
    }
}
