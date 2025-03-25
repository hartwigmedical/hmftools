package com.hartwig.hmftools.sage.tinc;

import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.LOGGER;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ALT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REF;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.GenotypeIds.fromVcfHeader;
import static com.hartwig.hmftools.common.variant.SageVcfTags.AVG_BASE_QUAL;
import static com.hartwig.hmftools.common.variant.SageVcfTags.MAP_QUAL_FACTOR;
import static com.hartwig.hmftools.common.variant.SageVcfTags.NEARBY_INDEL_FLAG;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_GERMLINE_ABQ_MIN;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_GERMLINE_DEPTH_HIGH;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_GERMLINE_DEPTH_LOW;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_MAX_FITTING_VARIANTS;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_MQF_LIMIT;
import static com.hartwig.hmftools.sage.tinc.TincConstants.TINC_RECOVERY_GERMLINE_AF_PROB;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
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

import org.apache.commons.math3.distribution.BinomialDistribution;

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

    public List<VariantData> variants() { return mVariants; }
    public List<VariantData> fittingVariants() { return mFittingVariants; }

    public void loadVariants()
    {
        if(mConfig.FitVariantsFile != null)
        {
            loadFitVariants();

            if(mConfig.RewriteVcf)
            {
                VcfFileReader vcfFileReader = new VcfFileReader(mConfig.InputVcf, true);

                for(VariantContext variantContext : vcfFileReader.iterator())
                {
                    mVariants.add(new VariantData(variantContext, mGenotypeIds));
                }

                SG_LOGGER.debug("loaded {} variants from VCF({})", mVariants.size(), mConfig.InputVcf);
            }

            return;
        }

        List<ChromosomeTask> chromosomeTasks = Lists.newArrayList();
        List<String> initialRefChromosomes = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

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

        LOGGER.debug("loading unfiltered VCF file({})", mConfig.InputVcf);

        final List<Callable> callableList = chromosomeTasks.stream().collect(Collectors.toList());

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
        {
            System.exit(1);
        }

        // calculate and then filter by GL_DP between 0.5x and 1.5x of mean (makes AD comparisons more meaningful with standardised DP)
        int[] filterCounts = new int[FilterReason.values().length];

        for(ChromosomeTask chrTask : chromosomeTasks)
        {
            mVariants.addAll(chrTask.variants());
            mFittingVariants.addAll(chrTask.fittingVariants());
            chrTask.addFilterCounts(filterCounts);
        }

        if(mFittingVariants.size() >= TINC_MAX_FITTING_VARIANTS)
        {
            downsampleFittingVariants();
        }

        filterByAverageReferenceDepth(filterCounts);

        String filterCountsStr = Arrays.stream(FilterReason.values()).map(x -> format("%s=%d", x.toString(), filterCounts[x.ordinal()]))
                .collect(Collectors.joining(","));

        LOGGER.debug("variants({}) fitCount({}) reasons: {}",
                mVariants.size(), mFittingVariants.size(), filterCountsStr);

        if(mConfig.WriteFitVariants)
            writeFitVariants();
    }

    private void downsampleFittingVariants()
    {
        int nthCount = (int)floor(mFittingVariants.size() / (double)TINC_MAX_FITTING_VARIANTS);

        if(nthCount < 2)
            return;

        int startCount = mFittingVariants.size();

        int index = 0;
        int counter = 0;
        while(index < mFittingVariants.size())
        {
            ++counter;

            if(counter < nthCount)
            {
                mFittingVariants.remove(index);
            }
            else
            {
                ++index;
                counter = 0;
            }
        }

        LOGGER.debug("down-sampled variants({} -> {})", startCount, mFittingVariants.size());
    }

    private void filterByAverageReferenceDepth(final int[] filterCounts)
    {
        long depthTotal = 0;

        for(VariantData variant : mFittingVariants)
        {
            depthTotal += variant.ReferenceDepth;
        }

        if(depthTotal == 0)
            return;

        int averageDepth = (int)round(depthTotal / (double)mFittingVariants.size());

        int index = 0;
        while(index < mFittingVariants.size())
        {
            VariantData variant = mFittingVariants.get(index);

            if(variant.ReferenceDepth < TINC_GERMLINE_DEPTH_LOW * averageDepth || variant.ReferenceDepth > TINC_GERMLINE_DEPTH_HIGH * averageDepth)
            {
                mFittingVariants.remove(index);
                ++filterCounts[FilterReason.REF_DEPTH.ordinal()];
            }
            else
            {
                ++index;
            }
        }

        if(filterCounts[FilterReason.REF_DEPTH.ordinal()] > 0)
        {
            SG_LOGGER.debug("filtered {} variants vs average ref depth({})",
                    filterCounts[FilterReason.REF_DEPTH.ordinal()], averageDepth);
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
        private final int[] mFilterCounts;

        public ChromosomeTask(final HumanChromosome chromosome)
        {
            mChromosome = chromosome;
            mChromosomeStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            mVariants = Lists.newArrayList();
            mFittingVariants = Lists.newArrayList();

            mFilterCounts = new int[FilterReason.values().length];
        }

        public List<VariantData> variants() { return mVariants; }
        public List<VariantData> fittingVariants() { return mFittingVariants; }

        public void addFilterCounts(final int[] filterCounts)
        {
            for(int i = 0; i < filterCounts.length; ++i)
            {
                filterCounts[i] += mFilterCounts[i];
            }
        }

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

            SG_LOGGER.trace("chr({}) starting variant annotation", mChromosome);

            for(VariantContext variantContext : vcfFileReader.regionIterator(chrRegion))
            {
                processVariant(variantContext);
                ++variantCount;

                if(variantCount > 0 && (variantCount % 100000) == 0)
                {
                    SG_LOGGER.debug("chr({}) processed {} variants", mChromosome, variantCount);
                }
            }

            LOGGER.debug("chr({}) loaded variants({}) unfiltered({})", mChromosome, mVariants.size(), mFittingVariants.size());

            mGnomadCache.removeCompleteChromosome(mChromosomeStr);
            mPonCache.removeCompleteChromosome(mChromosomeStr);

            return (long) 0;
        }

        private void processVariant(final VariantContext variantContext)
        {
            VariantData variant = new VariantData(variantContext, mGenotypeIds);
            mVariants.add(variant);

            setPonData(variant);

            FilterReason filterReason = checkFilters(variant);

            if(filterReason == FilterReason.NONE)
                mFittingVariants.add(variant);
            else
                ++mFilterCounts[filterReason.ordinal()];
        }

        @VisibleForTesting
        public static FilterReason checkFilters(final VariantData variant)
        {
            // annotate and evaluate for use as a fitting variant
            FilterReason filterReason = FilterReason.NONE;

            if(variant.isIndel())
            {
                filterReason = FilterReason.INDEL;
            }
            else if(!variant.isPassing() && !variant.OnlyGermlineFiltered)
            {
                filterReason = FilterReason.FILTERED;
            }
            else
            {
                // GL_ABQ >= 30 (removes low qual GL_AD variants from distorting the fit)
                int referenceABQ = Integer.parseInt(variant.RefGenotype.getAnyAttribute(AVG_BASE_QUAL).toString().split(CSV_DELIM)[1]);

                if(variant.ReferenceAltFrags > 0 && referenceABQ < TINC_GERMLINE_ABQ_MIN)
                    filterReason = FilterReason.LOW_ABQ;
            }

            if(filterReason == FilterReason.NONE)
            {
                double mqf = variant.Context.getAttributeAsDouble(MAP_QUAL_FACTOR, 0);

                if(mqf < TINC_MQF_LIMIT)
                    filterReason = FilterReason.MQF;
                else if(variant.Context.hasAttribute(NEARBY_INDEL_FLAG))
                    filterReason = FilterReason.NEAR_INDEL;
            }

            // check the PON for variants not otherwise filtered
            if(filterReason == FilterReason.NONE)
            {
                if(variant.gnomadFrequency() != null)
                {
                    filterReason = FilterReason.GNOMAD;
                }
                else if(variant.ponFiltered())
                {
                    filterReason = FilterReason.PON;
                }
            }

            if(filterReason == FilterReason.NONE && variant.ReferenceAltFrags > 0)
            {
                BinomialDistribution distribution = new BinomialDistribution(variant.ReferenceDepth, 0.5);

                if(distribution.cumulativeProbability(variant.ReferenceAltFrags) >= TINC_RECOVERY_GERMLINE_AF_PROB)
                {
                    filterReason = FilterReason.GERMLINE_AF;
                }
            }

            return filterReason;
        }

        private void setPonData(final VariantData variant)
        {
            if(mGnomadChrCache != null)
            {
                Double gnomadFreq = mGnomadChrCache.getFrequency(variant.isMnv(), variant.Ref, variant.Alt, variant.Position);

                if(gnomadFreq != null)
                {
                    variant.setGnomadFrequency(gnomadFreq);
                    variant.setPonFiltered();
                }
            }

            if(mPonChrCache != null)
            {
                PonVariantData ponData = mPonChrCache.getPonData(variant.Position, variant.Ref, variant.Alt);

                if(ponData != null)
                {
                    variant.setPonFrequency(ponData.Samples, ponData.MaxSampleReads, ponData.meanReadCount());

                    if(mPonCache.filterOnTierCriteria(variant.tier(), ponData.Samples, ponData.MaxSampleReads))
                        variant.setPonFiltered();
                }
            }
        }
    }

    private static final String FLD_REF_DEPTH = "RefDepth";
    private static final String FLD_TUMOR_DEPTH = "TumorDepth";
    private static final String FLD_REF_FRAGS = "RefFrags";
    private static final String FLD_TUMOR_FRAGS = "TumorFrags";

    private void writeFitVariants()
    {
        try
        {
            String filename = mConfig.InputVcf.replaceAll(".vcf.gz", ".fit_variant.tsv");
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME);
            sj.add(FLD_POSITION);
            sj.add(FLD_REF);
            sj.add(FLD_ALT);
            sj.add(FLD_REF_DEPTH);
            sj.add(FLD_REF_FRAGS);
            sj.add(FLD_TUMOR_DEPTH);
            sj.add(FLD_TUMOR_FRAGS);

            writer.write(sj.toString());
            writer.newLine();

            for(VariantData variant : mFittingVariants)
            {
                sj = new StringJoiner(TSV_DELIM);
                sj.add(variant.Chromosome);
                sj.add(String.valueOf(variant.Position));
                sj.add(variant.Ref);
                sj.add(variant.Alt);
                sj.add(String.valueOf(variant.ReferenceDepth));
                sj.add(String.valueOf(variant.ReferenceAltFrags));
                sj.add(String.valueOf(variant.TumorDepth));
                sj.add(String.valueOf(variant.TumorAltFrags));

                writer.write(sj.toString());
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write fit variants: {}", e.toString());
        }
    }

    private void loadFitVariants()
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.FitVariantsFile));

            String line = fileReader.readLine();
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, TSV_DELIM);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posIndex = fieldsIndexMap.get(FLD_POSITION);
            int refIndex = fieldsIndexMap.get(FLD_REF);
            int altIndex = fieldsIndexMap.get(FLD_ALT);
            int refDepthIndex = fieldsIndexMap.get(FLD_REF_DEPTH);
            int refFragsIndex = fieldsIndexMap.get(FLD_REF_FRAGS);
            int tumDepthIndex = fieldsIndexMap.get(FLD_TUMOR_DEPTH);
            int tumFragsIndex = fieldsIndexMap.get(FLD_TUMOR_FRAGS);

            while((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(TSV_DELIM, -1);

                VariantData variant = new VariantData(
                        values[chrIndex], Integer.parseInt(values[posIndex]), values[refIndex], values[altIndex],
                        Integer.parseInt(values[refDepthIndex]), Integer.parseInt(values[refFragsIndex]),
                        Integer.parseInt(values[tumDepthIndex]), Integer.parseInt(values[tumFragsIndex]));

                mFittingVariants.add(variant);
            }

            SG_LOGGER.info("loaded {} fit variants from file: {}", mFittingVariants.size(), mConfig.FitVariantsFile);
        }
        catch(IOException exception)
        {
            SG_LOGGER.error("failed to read fit variants file({})", mConfig.FitVariantsFile, exception.toString());
            System.exit(1);
        }
    }
}
