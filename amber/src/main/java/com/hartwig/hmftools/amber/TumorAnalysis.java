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

    public ListMultimap<Chromosome, TumorBAF> getBafs()
    {
        return mBafs;
    }

    public ListMultimap<Chromosome, TumorContamination> getContamination()
    {
        return mContamination;
    }

    public TumorAnalysis(
            final AmberConfig config, SamReaderFactory readerFactory,
            final ListMultimap<Chromosome, PositionEvidence> germlineHetLoci,
            final ListMultimap<Chromosome, PositionEvidence> germlineHomLoci)
            throws InterruptedException
    {
        mConfig = config;

        tumorBAFAndContamination(readerFactory, germlineHetLoci, germlineHomLoci);
    }

    private void tumorBAFAndContamination(final SamReaderFactory readerFactory,
            final ListMultimap<Chromosome, PositionEvidence> germlineHetLoci,
            final ListMultimap<Chromosome, PositionEvidence> germlineHomLoci) throws InterruptedException
    {
        AMB_LOGGER.info("processing tumor germline heterozygous({}) and homozygous({}) sites",
                germlineHetLoci.values().size(), germlineHomLoci.size());

        Map<Chromosome, List<PositionEvidence>> chrPositionEvidence = Maps.newHashMap();
        Map<PositionEvidence, PositionEvidence> contaminationBafMap = Maps.newHashMap();
        List<TumorBAF> tumorBAFs = Lists.newArrayList();

        for(Map.Entry<Chromosome, PositionEvidence> entry : germlineHetLoci.entries())
        {
            List<PositionEvidence> positions = chrPositionEvidence.get(entry.getKey());

            if(positions == null)
            {
                positions = Lists.newArrayList();
                chrPositionEvidence.put(entry.getKey(), positions);
            }

            PositionEvidence normal = entry.getValue();

            TumorBAF tumorBAF = TumorBAF.fromNormal(normal);
            tumorBAFs.add(tumorBAF);

            positions.add(tumorBAF.TumorEvidence);
        }

        for(Map.Entry<Chromosome, PositionEvidence> entry : germlineHomLoci.entries())
        {
            List<PositionEvidence> positions = chrPositionEvidence.get(entry.getKey());

            if(positions == null)
            {
                positions = Lists.newArrayList();
                chrPositionEvidence.put(entry.getKey(), positions);
            }

            PositionEvidence normal = entry.getValue();
            PositionEvidence tumor = PositionEvidence.copy(normal);

            positions.add(tumor);
            contaminationBafMap.put(normal, tumor);
        }

        // ensure positions are sorted after the merge
        for(List<PositionEvidence> positions : chrPositionEvidence.values())
        {
            Collections.sort(positions);
        }

        BamEvidenceReader bamEvidenceReader = new BamEvidenceReader(mConfig);
        bamEvidenceReader.processBam(mConfig.TumorBam, readerFactory, chrPositionEvidence);

        mBafs = ArrayListMultimap.create();

        tumorBAFs.stream()
                .filter(x -> x.TumorEvidence.IndelCount == 0)
                .forEach(x -> mBafs.put(HumanChromosome.fromString(x.chromosome()), x));

        mContamination = ArrayListMultimap.create();

        for(Map.Entry<PositionEvidence, PositionEvidence> entry : contaminationBafMap.entrySet())
        {
            PositionEvidence normal = entry.getKey();
            PositionEvidence tumor = entry.getValue();

            if(tumor.AltSupport > 0)
            {
                mContamination.put(
                        HumanChromosome.fromString(normal.chromosome()),
                        new TumorContamination(normal.Chromosome, normal.Position, normal.toBaseDepthData(), tumor.toBaseDepthData()));
            }
        }
        System.out.println(mContamination.size());
    }
}
