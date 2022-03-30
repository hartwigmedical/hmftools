package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthFactory;
import com.hartwig.hmftools.common.amber.ModifiableBaseDepth;
import com.hartwig.hmftools.common.amber.ModifiableTumorBAF;
import com.hartwig.hmftools.common.amber.TumorBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;

public class AmberTumor
{
    private final AmberConfig mConfig;
    private ListMultimap<Chromosome, TumorBAF> mBafs;
    private ListMultimap<Chromosome, TumorContamination> mContamination;

    public ListMultimap<Chromosome, TumorBAF> getBafs() { return mBafs; }
    public ListMultimap<Chromosome, TumorContamination> getContamination() { return mContamination; }

    AmberTumor(final AmberConfig config, SamReaderFactory readerFactory,
            final ListMultimap<Chromosome, BaseDepth> germlineHetLoci,
            final ListMultimap<Chromosome, BaseDepth> germlineHomLoci)
            throws InterruptedException
    {
        mConfig = config;

        tumorBAFAndContamination(readerFactory, germlineHetLoci, germlineHomLoci);
    }

    // we process them together
    private void tumorBAFAndContamination(final SamReaderFactory readerFactory,
            final ListMultimap<Chromosome, BaseDepth> germlineHetLoci, final ListMultimap<Chromosome, BaseDepth> germlineHomLoci) throws InterruptedException
    {
        AMB_LOGGER.info("Processing {} germline heterozygous loci in tumor bam {}", germlineHetLoci.values().size(), mConfig.TumorBamPath);
        AMB_LOGGER.info("Processing {} germline homozygous loci in tumor bam {} for contamination", germlineHomLoci.size(), mConfig.TumorBamPath);

        final List<ModifiableTumorBAF> tumorBAFs = germlineHetLoci.values().stream().sorted().map(TumorBAFFactory::create).collect(Collectors.toList());
        final Map<BaseDepth, ModifiableBaseDepth> contaminationBafMap = germlineHomLoci.values().stream().collect(
                Collectors.toMap(x -> x, BaseDepthFactory::create));

        var tumorBafFactory = new TumorBAFFactory(mConfig.MinBaseQuality);
        var contaminationBafFactory = new BaseDepthFactory(mConfig.MinBaseQuality);

        // merge both into one list
        final List<GenomePosition> mergedLoci = new ArrayList<>(tumorBAFs);
        mergedLoci.addAll(contaminationBafMap.values());
        mergedLoci.sort(null);

        BiConsumer<GenomePosition, SAMRecord> lociBamRecordHander = (GenomePosition genomePosition, SAMRecord samRecord) ->
        {
            if (genomePosition instanceof ModifiableTumorBAF)
                tumorBafFactory.addEvidence((ModifiableTumorBAF)genomePosition, samRecord);
            else if (genomePosition instanceof ModifiableBaseDepth)
                contaminationBafFactory.addEvidence((ModifiableBaseDepth)genomePosition, samRecord);
        };

        AsyncBamLociReader.processBam(mConfig.TumorBamPath, readerFactory, mergedLoci,
                lociBamRecordHander, mConfig.ThreadCount, mConfig.MinMappingQuality);

        mBafs = ArrayListMultimap.create();

        tumorBAFs.stream().filter(x -> x.tumorIndelCount() == 0).forEach(x ->
            mBafs.put(HumanChromosome.fromString(x.chromosome()), x));

        mContamination = ArrayListMultimap.create();

        contaminationBafMap.forEach((normal, tumor) ->
        {
            if (tumor.altSupport() != 0)
            {
                mContamination.put(HumanChromosome.fromString(normal.chromosome()),
                        ImmutableTumorContamination.builder().from(normal).normal(normal).tumor(tumor).build());
            }
        });
    }
}
