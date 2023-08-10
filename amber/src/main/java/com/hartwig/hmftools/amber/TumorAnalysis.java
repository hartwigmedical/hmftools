package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import htsjdk.samtools.SamReaderFactory;

public class TumorAnalysis
{
    private final AmberConfig mConfig;
    private ListMultimap<Chromosome, TumorBAF> mBafs;
    private ListMultimap<Chromosome, TumorContamination> mContamination;

    public ListMultimap<Chromosome, TumorBAF> getBafs() { return mBafs; }
    public ListMultimap<Chromosome, TumorContamination> getContamination() { return mContamination; }

    public TumorAnalysis(
            final AmberConfig config, SamReaderFactory readerFactory,
            final ListMultimap<Chromosome, BaseDepth> germlineHetLoci,
            final ListMultimap<Chromosome, BaseDepth> germlineHomLoci)
            throws InterruptedException
    {
        mConfig = config;

        tumorBAFAndContamination(readerFactory, germlineHetLoci, germlineHomLoci);
    }

    // we process them together
    private void tumorBAFAndContamination(final SamReaderFactory readerFactory,
            final ListMultimap<Chromosome,BaseDepth> germlineHetLoci, final ListMultimap<Chromosome, BaseDepth> germlineHomLoci) throws InterruptedException
    {
        AMB_LOGGER.info("processing {} germline heterozygous sites", germlineHetLoci.values().size());
        AMB_LOGGER.info("processing {} germline homozygous sites", germlineHomLoci.size());

        /*
        final List<TumorBAF> tumorBAFs = germlineHetLoci.values().stream().sorted().map(TumorBAFFactory::create).collect(Collectors.toList());
        final Map<BaseDepth,BaseDepth> contaminationBafMap = germlineHomLoci.values().stream().collect(
                Collectors.toMap(x -> x, BaseDepthFactory::copyBaseDepth));
        */

        Map<Chromosome,List<BaseDepth>> chrBaseDepthMap = Maps.newHashMap();
        Map<BaseDepth,BaseDepth> contaminationBafMap = Maps.newHashMap();
        List<TumorBAF> tumorBAFs = Lists.newArrayList();

        for(Map.Entry<Chromosome,BaseDepth> entry : germlineHetLoci.entries())
        {
            List<BaseDepth> positions = chrBaseDepthMap.get(entry.getKey());

            if(positions == null)
            {
                positions = Lists.newArrayList();
                chrBaseDepthMap.put(entry.getKey(), positions);
            }

            BaseDepth normal = entry.getValue();

            TumorBAF tumorBAF = TumorBAFFactory.create(normal);

            positions.add(tumorBAF.TumorEvidence);
        }

        for(Map.Entry<Chromosome,BaseDepth> entry : germlineHomLoci.entries())
        {
            List<BaseDepth> positions = chrBaseDepthMap.get(entry.getKey());

            if(positions == null)
            {
                positions = Lists.newArrayList();
                chrBaseDepthMap.put(entry.getKey(), positions);
            }

            BaseDepth normal = entry.getValue();
            BaseDepth tumor = BaseDepth.copy(normal);

            positions.add(tumor);
            contaminationBafMap.put(normal, tumor);
        }

        // ensure positions are sorted after the merge
        for(List<BaseDepth> positions : chrBaseDepthMap.values())
        {
            Collections.sort(positions);
        }

        BamEvidenceReader bamEvidenceReader = new BamEvidenceReader(mConfig);
        bamEvidenceReader.processBam(mConfig.TumorBam, readerFactory, chrBaseDepthMap);

        // OLD CODE, to be deprecated

        /*
            TumorBAFFactory tumorBafFactory = new TumorBAFFactory(mConfig.MinBaseQuality);
            BaseDepthFactory contaminationBafFactory = new BaseDepthFactory(mConfig.MinBaseQuality);

            // merge both into one list
            final List<GenomePosition> mergedLoci = new ArrayList<>(tumorBAFs);
            mergedLoci.addAll(contaminationBafMap.values());
            mergedLoci.sort(null);

            BiConsumer<GenomePosition, SAMRecord> lociBamRecordHander = (GenomePosition genomePosition, SAMRecord samRecord) ->
            {
                if (genomePosition instanceof TumorBAF)
                    tumorBafFactory.addEvidence((TumorBAF)genomePosition, samRecord);
                else if (genomePosition instanceof BaseDepth)
                    contaminationBafFactory.addEvidence((BaseDepth)genomePosition, samRecord);
            };

            AsyncBamLociReader.processBam(
                    mConfig.TumorBam, readerFactory, mergedLoci, lociBamRecordHander, mConfig.Threads, mConfig.MinMappingQuality);
        */

        mBafs = ArrayListMultimap.create();

        tumorBAFs.stream().filter(x -> x.TumorIndelCount == 0).forEach(x -> mBafs.put(HumanChromosome.fromString(x.chromosome()), x));

        mContamination = ArrayListMultimap.create();

        contaminationBafMap.forEach((normal, tumor) ->
        {
            if(tumor.AltSupport != 0)
            {
                TumorContamination tumorContamination = ImmutableTumorContamination.builder()
                        .chromosome(normal.chromosome())
                        .position(normal.position())
                        .normal(normal.toBaseDepthData())
                        .tumor(tumor.toBaseDepthData())
                        .build();

                mContamination.put(HumanChromosome.fromString(normal.chromosome()), tumorContamination);
            }
        });
    }
}
