package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.paneldesign.DataWriter.CANDIDATE_FILE_EXTENSION;
import static com.hartwig.hmftools.geneutils.paneldesign.DataWriter.GENE_REGION_FILE_EXTENSION;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.GENE_CANDIDATE_REGION_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.GENE_FLANKING_DISTANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.GENE_LONG_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.GENE_MAX_CANDIDATE_PROBES;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.GENE_MAX_EXONS_TO_ADD_INTRON;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.GENE_MIN_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.MIN_PROBE_QUALITY_SCORE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.PROBE_GC_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.PROBE_GC_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.ProbeCandidate.createProbeCandidate;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

public class GeneProbesGenerator
{
    private static class GeneNameTranscriptId
    {
        public final String GeneName;
        public final String TranscriptName;

        public GeneNameTranscriptId(final String geneName, final String transcriptName)
        {
            GeneName = geneName;
            TranscriptName = transcriptName;
        }

        public boolean useCanonical() { return TranscriptName.isEmpty(); }
        public String toString() { return format("%s:%s", GeneName, useCanonical() ? "canonical" : TranscriptName); }
    }

    private final PanelConfig mConfig;
    private final PanelCache mPanelCache;
    private final ProbeQualityProfile mProbeQualityProfile;

    private final EnsemblDataCache mEnsemblDataCache;

    private final List<GeneNameTranscriptId> mGeneNameTranscriptIds;

    public GeneProbesGenerator(final PanelConfig config, final PanelCache panelCache, final ProbeQualityProfile probeQualityProfile)
    {
        mConfig = config;
        mPanelCache = panelCache;
        mProbeQualityProfile = probeQualityProfile;

        mGeneNameTranscriptIds = new ArrayList<>();

        mEnsemblDataCache = new EnsemblDataCache(mConfig.EnsemblDir, mConfig.RefGenVersion);
        mEnsemblDataCache.setRequiredData(true, false, false, false);
        mEnsemblDataCache.load(false);

        try(DelimFileReader reader = new DelimFileReader(mConfig.GeneTranscriptFile))
        {
            reader.stream().forEach(row -> mGeneNameTranscriptIds.add(
                    new GeneNameTranscriptId(row.get(FLD_GENE_NAME), row.getOrNull(FLD_TRANS_NAME))));
        }
    }

    public void run()
    {
        List<TargetedGene> targetedGenes = new ArrayList<>();

        // from the gene, find all the gene regions
        for(GeneNameTranscriptId geneTranscript : mGeneNameTranscriptIds)
        {
            String geneName = geneTranscript.GeneName;

            GU_LOGGER.debug("processing gene-transcript({})", geneTranscript);

            // get all the gene regions from the ensembl
            GeneData geneData = mEnsemblDataCache.getGeneDataByName(geneName);

            if(geneData == null)
            {
                GU_LOGGER.error("transcript data for gene: {} not found", geneName);
                throw new RuntimeException(format("transcript data for gene: %s not found", geneName));
            }

            // skip if the gene is already covered by a custom region
            if(mPanelCache.overlapsExisting(geneData.Chromosome, geneData.GeneStart, geneData.GeneEnd))
            {
                GU_LOGGER.warn("gene({}) skipped since overlaps with existing panel region", geneName);
                continue;
            }

            TranscriptData transcriptData = geneTranscript.useCanonical() ?
                    mEnsemblDataCache.getCanonicalTranscriptData(geneData.GeneId) :
                    mEnsemblDataCache.getTranscriptData(geneData.GeneId, geneTranscript.TranscriptName);

            if(transcriptData == null)
            {
                GU_LOGGER.error("gene({}) transcript({}) not found", geneName, geneTranscript.TranscriptName);
                System.exit(1);
            }

            if(transcriptData.nonCoding())
            {
                // should add as a custom region instead
                GU_LOGGER.warn("gene({}) transcript({}) non-coding skipped", geneName, geneTranscript.TranscriptName);
                continue;
            }

            TargetedGene targetedGene = new TargetedGene(geneData, transcriptData);
            targetedGenes.add(targetedGene);
        }

        for(TargetedGene targetedGene : targetedGenes)
        {
            populateTargetedGeneRegions(targetedGene);
            populateCandidateProbes(targetedGene);
        }

        computeCandidateProbeQualityScores(targetedGenes);

        // now choose probes for each region
        selectProbeCandidates(targetedGenes);

        // add regions & probes to the cache
        List<TargetedGeneRegion> targetedGeneRegions = Lists.newArrayList();

        for(TargetedGene targetedGene : targetedGenes)
        {
            for(TargetedGeneRegion targetedGeneRegion : targetedGene.getRegions())
            {
                ChrBaseRegion region = new ChrBaseRegion(
                        targetedGene.getGeneData().Chromosome, targetedGeneRegion.getStart(), targetedGeneRegion.getEnd());

                PanelRegion panelRegion;

                String sourceInfo = format("%s:%s", targetedGene.getGeneData().GeneName, targetedGeneRegion.getType().toString());

                if(targetedGeneRegion.useWholeRegion())
                {
                    panelRegion = new PanelRegion(region, RegionType.GENE, sourceInfo);
                }
                else
                {
                    ProbeCandidate probeCandidate = targetedGeneRegion.getSelectedProbe();

                    if(probeCandidate == null)
                        continue;

                    panelRegion = new PanelRegion(
                            region, RegionType.GENE, sourceInfo,
                            probeCandidate.getSequence(), probeCandidate.getGcContent(), probeCandidate.getQualityScore().get());
                }

                mPanelCache.addRegion(panelRegion);

                targetedGeneRegions.add(targetedGeneRegion);
            }
        }

        // write the outputs
        DataWriter.writeCandidates(mConfig.formOutputFilename(CANDIDATE_FILE_EXTENSION), targetedGeneRegions);
        DataWriter.writeTargertedGeneRegions(mConfig.formOutputFilename(GENE_REGION_FILE_EXTENSION), targetedGeneRegions);
    }

    /*
    | Region Type | Gene                  | Design                                                                                         |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    | Coding      | All                   | - Cover the full coding region of each exon                                                    |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    | UTR         | All                   | - Add a 120nt probe centred on each non-coding exon                                            |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    | Upstream/   | All                   | - Create 8 candidate 120nt probes from 1kb to 2kb bases upstream of gene.                      |
    | Downstream  |                       | - Filter for  0.35<GC<0.55 AND SUM_BLASTN<1000                                                 |
    |             |                       | - Choose 1 probe per region with lowest SUM_BLASTN score.                                      |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    | Intronic    | Genes with < 20 exons | - For each intron of >5kb, create 8 candidate120nt probes from 1kb-2kb from EACH flanking exon |
    |             |                       | - For each intron of 3kb-5kb, create 8 candidate 120nt probes centred on intron                |
    |             |                       | - Filter for 0.35<GC<0.55 AND SUM_BLASTN<1000                                                  |
    |             |                       | - Choose 1 probe  per region with lowest SUM_BLASTN score.                                     |
    |-------------|-----------------------|------------------------------------------------------------------------------------------------|
    */
    static void populateTargetedGeneRegions(final TargetedGene targetedGene)
    {
        // first we create a region 1-2kb upstream of the gene
        TranscriptData transcript = targetedGene.getTranscriptData();

        targetedGene.addRegion(targetedGene.getGeneData().forwardStrand() ? TargetedGeneRegion.Type.UP_STREAM : TargetedGeneRegion.Type.DOWN_STREAM,
                transcript.TransStart - GENE_FLANKING_DISTANCE - GENE_CANDIDATE_REGION_SIZE,
                transcript.TransStart - GENE_FLANKING_DISTANCE - 1);

        int lastExonEnd = -1;

        for(ExonData exonData : transcript.exons())
        {
            if(lastExonEnd != -1 && targetedGene.getTranscriptData().exons().size() <= GENE_MAX_EXONS_TO_ADD_INTRON)
            {
                int intronLength = exonData.Start - lastExonEnd;

                if(intronLength > GENE_MIN_INTRON_LENGTH)
                {
                    if(intronLength > GENE_LONG_INTRON_LENGTH)
                    {
                        // for long intron, we want to split into two 1k regions each flank the exons
                        targetedGene.addRegion(TargetedGeneRegion.Type.INTRONIC_LONG,
                                lastExonEnd + 1 + GENE_FLANKING_DISTANCE,
                                lastExonEnd + GENE_FLANKING_DISTANCE + GENE_CANDIDATE_REGION_SIZE);

                        targetedGene.addRegion(TargetedGeneRegion.Type.INTRONIC_LONG,
                                exonData.Start - GENE_FLANKING_DISTANCE - GENE_CANDIDATE_REGION_SIZE,
                                exonData.Start - GENE_FLANKING_DISTANCE - 1);
                    }
                    else
                    {
                        // for smaller intron we add a region in the middle 1k
                        int intronMid = (lastExonEnd + 1 + exonData.Start) / 2;
                        targetedGene.addRegion(TargetedGeneRegion.Type.INTRONIC_SHORT,
                                intronMid - GENE_CANDIDATE_REGION_SIZE / 2,
                                intronMid + GENE_CANDIDATE_REGION_SIZE / 2 - 1);
                    }
                }
            }

            boolean isCoding = positionsOverlap(exonData.Start, exonData.End, transcript.CodingStart, transcript.CodingEnd);

            if(isCoding)
            {
                // limit the region's bases to coding
                targetedGene.addRegion(TargetedGeneRegion.Type.CODING,
                        Math.max(exonData.Start, transcript.CodingStart),
                        Math.min(exonData.End, transcript.CodingEnd));
            }
            else
            {
                // just one probe at the middle, so just the 120 bases
                int exonMid = (exonData.Start + exonData.End + 1) / 2;
                targetedGene.addRegion(TargetedGeneRegion.Type.UTR, exonMid - PROBE_LENGTH / 2, exonMid + PROBE_LENGTH / 2 - 1);
            }

            lastExonEnd = exonData.End;
        }

        targetedGene.addRegion(targetedGene.getGeneData().forwardStrand() ? TargetedGeneRegion.Type.DOWN_STREAM : TargetedGeneRegion.Type.UP_STREAM,
                transcript.TransEnd + 1 + GENE_FLANKING_DISTANCE,
                transcript.TransEnd + GENE_FLANKING_DISTANCE + GENE_CANDIDATE_REGION_SIZE);
    }

    // from the targeted gene regions, we find the probes
    void populateCandidateProbes(TargetedGene targetedGene)
    {
        for(TargetedGeneRegion region : targetedGene.getRegions())
        {
            populateCandidateProbes(region, mConfig.RefGenome);
        }
    }

    static void populateCandidateProbes(final TargetedGeneRegion region, final RefGenomeInterface refGenomeInterface)
    {
        switch(region.getType())
        {
            case CODING:
            case UTR:
                // use whole region, so do not add probe candidate
                break;
            case UP_STREAM:
            case DOWN_STREAM:
            case INTRONIC_LONG:
            case INTRONIC_SHORT:

                // create set of candidate probes within the region
                for(int i = 0; i < GENE_MAX_CANDIDATE_PROBES; ++i)
                {
                    int start = region.getStart() + i * PROBE_LENGTH;
                    int end = start + PROBE_LENGTH - 1; // end is inclusive

                    ProbeCandidate probe = createProbeCandidate(
                            new ChrBaseRegion(region.getGene().getGeneData().Chromosome, start, end), refGenomeInterface);

                    region.getProbeCandidates().add(probe);
                }
                break;
        }
    }

    public void computeCandidateProbeQualityScores(final List<TargetedGene> targetedGeneList)
    {
        List<ProbeCandidate> probeCandidates = Lists.newArrayList();

        for(TargetedGene targetedGene : targetedGeneList)
        {
            for(TargetedGeneRegion targetedGeneRegion : targetedGene.getRegions())
            {
                for(ProbeCandidate probeCandidate : targetedGeneRegion.getProbeCandidates())
                {
                    if(checkGcContent(probeCandidate))
                    {
                        probeCandidates.add(probeCandidate);
                    }
                }
            }
        }

        if(probeCandidates.isEmpty())
            return;

        GU_LOGGER.info("Computing quality scores of gene probes");
        List<Optional<Double>> qualityScores = probeCandidates.stream()
                .map(probe -> mProbeQualityProfile.computeQualityScore(probe.region())).toList();

        for(int i = 0; i < probeCandidates.size(); ++i)
        {
            ProbeCandidate probeCandidate = probeCandidates.get(i);
            qualityScores.get(i).ifPresent(probeCandidate::setQualityScore);
        }
    }

    public void selectProbeCandidates(Collection<TargetedGene> targetedGeneList)
    {
        // select top-scoring, valid probe from the set of candidates
        for(TargetedGene targetedGene : targetedGeneList)
        {
            for(TargetedGeneRegion targetedGeneRegion : targetedGene.getRegions())
            {
                ProbeCandidate selectedProbe = null;

                for(ProbeCandidate probeCandidate : targetedGeneRegion.getProbeCandidates())
                {
                    checkSetProbeCandidateFilter(probeCandidate);

                    if(!probeCandidate.passFilter())
                    {
                        continue;
                    }

                    if(probeCandidate.getQualityScore().isPresent() &&
                       (selectedProbe == null || probeCandidate.getQualityScore().get() > selectedProbe.getQualityScore().get()))
                    {
                        selectedProbe = probeCandidate;
                    }
                }

                if(selectedProbe != null)
                {
                    targetedGeneRegion.setSelectedProbe(selectedProbe);

                    GU_LOGGER.trace("region: type({}) {}:{}-{} selected probe: {}-{}", targetedGeneRegion.getType(),
                            targetedGeneRegion.getChromosome(), targetedGeneRegion.getStart(), targetedGeneRegion.getEnd(),
                            selectedProbe.getStart(), selectedProbe.getEnd());
                }
            }
        }
    }

    public static boolean checkGcContent(ProbeCandidate probeCandidate)
    {
        return probeCandidate.getGcContent() > PROBE_GC_MIN && probeCandidate.getGcContent() < PROBE_GC_MAX;
    }

    public static void checkSetProbeCandidateFilter(ProbeCandidate probeCandidate)
    {
        if(probeCandidate.getGcContent() <= PROBE_GC_MIN)
        {
            probeCandidate.setFilterReason(format("gc <= %.2f", PROBE_GC_MIN));
        }
        else if(probeCandidate.getGcContent() >= PROBE_GC_MAX)
        {
            probeCandidate.setFilterReason(format("gc >= %.2f", PROBE_GC_MAX));
        }
        else if(probeCandidate.getQualityScore().orElse(0d) < MIN_PROBE_QUALITY_SCORE)
        {
            probeCandidate.setFilterReason(format("quality < %g", MIN_PROBE_QUALITY_SCORE));
        }
    }
}
