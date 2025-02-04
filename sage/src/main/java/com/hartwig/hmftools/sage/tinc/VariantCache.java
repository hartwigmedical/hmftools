package com.hartwig.hmftools.sage.tinc;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_BASE_QUAL;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.tinc.TincCalculator.TINC_GERMLINE_ABQ_MIN;
import static com.hartwig.hmftools.sage.tinc.TincCalculator.TINC_GERMLINE_DEPTH_HIGH;
import static com.hartwig.hmftools.sage.tinc.TincCalculator.TINC_GERMLINE_DEPTH_LOW;
import static com.hartwig.hmftools.sage.tinc.TincCalculator.TINC_MAX_FITTING_VARIANTS;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.pon.GnomadCache;
import com.hartwig.hmftools.common.variant.pon.GnomadChrCache;
import com.hartwig.hmftools.common.variant.pon.PonCache;
import com.hartwig.hmftools.common.variant.pon.PonChrCache;
import com.hartwig.hmftools.common.variant.pon.PonVariantData;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VariantCache
{
    private final TincConfig mConfig;
    private final GnomadCache mGnomadCache;
    private final PonCache mPonCache;

    private final List<VariantData> mVariants;
    private final List<VariantData> mFittingVariants;

    private final GenotypeIds mGenotypeIds;

    private static final List<String> PERMITTED_FILTERS = Lists.newArrayList(
            "maxGermlineRelQual", "maxGermlineVAF", "maxGermlineRelQual", "maxGermlineVAF", PASS);

    public VariantCache(final TincConfig config)
    {
        mConfig = config;

        mVariants = Lists.newArrayList();
        mFittingVariants = Lists.newArrayList();

        mGnomadCache = new GnomadCache(mConfig.RefGenVersion, mConfig.GnomadFile, mConfig.GnomadDirectory);
        mPonCache = new PonCache(mConfig.PonFilename, true);

        if(Files.exists(Paths.get(mConfig.InputVcf)))
        {
            VcfFileReader vcfFileReader = new VcfFileReader(mConfig.InputVcf);
            VCFHeader vcfHeader = vcfFileReader.vcfHeader();
            mGenotypeIds = fromVcfHeader(vcfHeader, mConfig.ReferenceId, mConfig.TumorId);
        }
        else
        {
            mGenotypeIds = null;
        }
    }

    public void loadVariants()
    {
        List<ChromosomeTask> chromosomeTasks = Lists.newArrayList();
        List<String> initialRefChromosomes = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            /*
            if(!mConfig.SpecificRegions.isEmpty() && mConfig.SpecificRegions.stream().noneMatch(x -> x.Chromosome.equals(chrStr)))
            {
                mVcfWriter.onChromosomeComplete(chromosome);
                continue;
            }
            */

            ChromosomeTask chromosomeTask = new ChromosomeTask(chromosome);

            chromosomeTasks.add(chromosomeTask);

            if(initialRefChromosomes.size() < max(mConfig.Threads, 1))
                initialRefChromosomes.add(chrStr);
        }

        for(String chromosome : initialRefChromosomes)
        {
            mPonCache.loadPonEntries(chromosome);
        }

        mGnomadCache.initialise(initialRefChromosomes);

        LOGGER.info("loading unfiltered VCF file({})", mConfig.InputVcf);

        final List<Callable> callableList = chromosomeTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
        {
            System.exit(1);
        }

        // calculate and then filter by GL_DP between 0.5x and 1.5x of mean (makes AD comparisons more meaningful with standardised DP)

        for(ChromosomeTask chrTask : chromosomeTasks)
        {
            mVariants.addAll(chrTask.variants());
            mFittingVariants.addAll(chrTask.fittingVariants());
        }

        LOGGER.info("variants({}) filtered({})", mVariants.size(), mFittingVariants.size());

        if(mFittingVariants.size() >= TINC_MAX_FITTING_VARIANTS)
        {
            downsampleFittingVariants();
        }

        filterByAverageReferenceDepth();
    }

    private void downsampleFittingVariants()
    {
        int nthCount = (int)floor(mFittingVariants.size() / (double)TINC_MAX_FITTING_VARIANTS);

        if(nthCount < 2)
            return;

        int index = 0;
        int counter = 0;
        while(index < mFittingVariants.size())
        {
            ++counter;

            if(counter < nthCount)
            {
                mFittingVariants.remove(index);
                counter = 0;
            }
            else
            {
                ++index;
            }
        }

        LOGGER.info("down-sampled variants({})", mFittingVariants.size());
    }

    private void filterByAverageReferenceDepth()
    {
        long depthTotal = 0;

        for(VariantData variant : mFittingVariants)
        {
            int refDepth = variant.Context.getGenotype(mGenotypeIds.ReferenceOrdinal).getDP();
            depthTotal += refDepth;
        }

        if(depthTotal == 0)
            return;

        int averageDepth = (int)round(depthTotal / (double)mFittingVariants.size());

        int filtered = 0;

        int index = 0;
        while(index < mFittingVariants.size())
        {
            VariantData variant = mFittingVariants.get(index);

            int refDepth = variant.Context.getGenotype(mGenotypeIds.ReferenceOrdinal).getDP();

            if(refDepth < TINC_GERMLINE_DEPTH_LOW * averageDepth || refDepth > TINC_GERMLINE_DEPTH_HIGH * averageDepth)
            {
                mFittingVariants.remove(index);
                ++filtered;
            }
            else
            {
                ++index;
            }
        }

        if(filtered > 0)
        {
            SG_LOGGER.debug("filtered {} variants vs average ref depth({})", filtered, averageDepth);
        }
    }

    public class ChromosomeTask implements Callable
    {
        private final HumanChromosome mChromosome;
        private final String mChromosomeStr;

        private final List<VariantData> mVariants;
        private final List<VariantData> mFittingVariants;

        private GnomadChrCache mGnomadChrCache;
        private PonChrCache mPonChrCache;

        public ChromosomeTask(final HumanChromosome chromosome)
        {
            mChromosome = chromosome;
            mChromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            mVariants = Lists.newArrayList();
            mFittingVariants = Lists.newArrayList();
        }

        public List<VariantData> variants() { return mVariants; }
        public List<VariantData> fittingVariants() { return mFittingVariants; }

        @Override
        public Long call()
        {
            int variantCount = 0;

            VcfFileReader vcfFileReader = new VcfFileReader(mConfig.InputVcf, true);

            if(!vcfFileReader.fileValid())
            {
                SG_LOGGER.error("invalid somatic VCF file({})", mConfig.InputVcf);
                System.exit(1);
            }

            mGnomadChrCache = mGnomadCache.getChromosomeCache(mChromosomeStr);
            mPonChrCache = mPonCache.getChromosomeCache(mChromosomeStr);

            RefGenomeCoordinates coordinates =
                    mConfig.RefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
            ChrBaseRegion chrRegion = new ChrBaseRegion(mChromosomeStr, 1, coordinates.Lengths.get(mChromosome));

            SG_LOGGER.debug("chr({}) starting variant annotation", mChromosome);

            for(VariantContext variantContext : vcfFileReader.regionIterator(chrRegion))
            {
                /*
            }
                if(!mConfig.SpecificRegions.isEmpty())
                {
                    if(mConfig.SpecificRegions.stream()
                            .noneMatch(x -> x.containsPosition(variantContext.getContig(), variantContext.getStart())))
                        continue;
                }
                */

                processVariant(variantContext);
                ++variantCount;

                if(variantCount > 0 && (variantCount % 100000) == 0)
                {
                    SG_LOGGER.debug("chr({}) processed {} variants", mChromosome, variantCount);
                }
            }

            return (long) 0;
        }

        private void processVariant(final VariantContext variantContext)
        {
            VariantData variant = new VariantData(variantContext);
            mVariants.add(variant);

            if(variant.isIndel())
                return;

            if(variant.Context.getFilters().stream().anyMatch(x -> !PERMITTED_FILTERS.contains(x)))
                return;

            // annotate and evaluate for use as a fitting variant
            Double gnomadFreq = mGnomadChrCache.getFrequency(variant.isMnv(), variant.Ref, variant.Alt, variant.Position);

            if(gnomadFreq != null)
                return;

            if(mPonChrCache != null)
            {
                PonVariantData ponData = mPonChrCache.getPonData(variant.Position, variant.Ref, variant.Alt);

                if(ponData != null && mPonCache.filterOnTierCriteria(variant.tier(), ponData.Samples, ponData.MaxSampleReads))
                    return;
            }

            // GL_ABQ >= 30 (removes low qual 1 GL_AD variants from distorting the fit).
            // TODO: Note that we currently don't annotate this field (defaults to 0) if GL_AD = 0, so we should change this
            Genotype refGenotype = variant.Context.getGenotype(mGenotypeIds.ReferenceOrdinal);
            int germlineABQ = Integer.parseInt(refGenotype.getAnyAttribute(AVG_BASE_QUAL).toString().split(CSV_DELIM)[1]);

            if(germlineABQ < TINC_GERMLINE_ABQ_MIN)
                return;

            mFittingVariants.add(variant);
        }
    }
}
