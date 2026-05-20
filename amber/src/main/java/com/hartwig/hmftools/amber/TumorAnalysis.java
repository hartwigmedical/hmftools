package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.amber.AmberConstants.CONTAMINATON_MIN_NORMAL_READ_DEPTH;
import static com.hartwig.hmftools.amber.AmberUtils.aboveQualFilter;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.amber.contamination.TumorContamination;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import htsjdk.samtools.SamReaderFactory;

public class TumorAnalysis
{
    private final AmberConfig mConfig;
    private final ListMultimap<Chromosome, TumorBAF> mBafs;
    private final List<TumorContamination> mContaminationSites;

    public ListMultimap<Chromosome, TumorBAF> chrBafMap()
    {
        return mBafs;
    }
    public List<TumorContamination> contaminationSites() { return mContaminationSites; }

    public TumorAnalysis(
            final AmberConfig config, SamReaderFactory readerFactory,
            final ListMultimap<Chromosome, PositionEvidence> germlineHetLoci,
            final ListMultimap<Chromosome, PositionEvidence> germlineHomLoci)
    {
        mConfig = config;
        mContaminationSites = Lists.newArrayList();
        mBafs = ArrayListMultimap.create();

        tumorBAFAndContamination(readerFactory, germlineHetLoci, germlineHomLoci);
    }

    private void tumorBAFAndContamination(
            final SamReaderFactory readerFactory, final ListMultimap<Chromosome, PositionEvidence> germlineHetLoci,
            final ListMultimap<Chromosome, PositionEvidence> germlineHomLoci)
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

        tumorBAFs.stream()
                .filter(x -> x.TumorEvidence.IndelCount == 0)
                .forEach(x -> mBafs.put(HumanChromosome.fromString(x.chromosome()), x));

        for(Map.Entry<PositionEvidence, PositionEvidence> entry : contaminationBafMap.entrySet())
        {
            PositionEvidence normal = entry.getKey();

            if(normal.RefSupport <= CONTAMINATON_MIN_NORMAL_READ_DEPTH)
                continue;

            PositionEvidence tumor = entry.getValue();

            if(tumor.AltSupport > 0 && aboveQualFilter(tumor))
            {
                mContaminationSites.add(
                        new TumorContamination(normal.Chromosome, normal.Position, normal.toBaseDepthData(), tumor.toBaseDepthData()));
            }
        }
    }
}
