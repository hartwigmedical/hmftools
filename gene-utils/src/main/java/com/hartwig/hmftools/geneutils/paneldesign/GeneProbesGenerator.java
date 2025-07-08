package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.paneldesign.DataWriter.CANDIDATE_FILE_EXTENSION;
import static com.hartwig.hmftools.geneutils.paneldesign.DataWriter.GENE_REGION_FILE_EXTENSION;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_CANDIDATE_REGION_SIZE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_FLANKING_DISTANCE;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_LONG_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MAX_CANDIDATE_PROBES;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MAX_EXONS_TO_ADD_INTRON;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.GENE_MIN_INTRON_LENGTH;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_QUALITY_SCORE_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_GC_MAX;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_GC_MIN;
import static com.hartwig.hmftools.geneutils.paneldesign.PanelBuilderConstants.PROBE_LENGTH;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;

public class GeneProbesGenerator
{
    private final PanelBuilderConfig mConfig;
    private final PanelCache mPanelCache;
    private final ProbeQualityProfile mProbeQualityProfile;

    private final EnsemblDataCache mEnsemblDataCache;

    private final List<GeneNameTranscriptId> mGeneNameTranscriptIds;

    public GeneProbesGenerator(final PanelBuilderConfig config, final PanelCache panelCache, final ProbeQualityProfile probeQualityProfile)
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

        for(Gene targetedGene : targetedGenes)
        {
            populateTargetedGeneRegions(targetedGene);
            populateCandidateProbes(targetedGene);
        }

        computeCandidateProbeQualityScores(targetedGenes);

        // now choose probes for each region
        selectProbeCandidates(targetedGenes);

        // add regions & probes to the cache
        List<GeneRegion> targetedGeneRegions = Lists.newArrayList();

        for(Gene targetedGene : targetedGenes)
        {
            for(GeneRegion targetedGeneRegion : targetedGene.getRegions())
            {
                ChrBaseRegion region = new ChrBaseRegion(
                        targetedGene.getGeneData().Chromosome, targetedGeneRegion.getStart(), targetedGeneRegion.getEnd());

                ProbeCandidate panelRegion;

                String sourceInfo = format("%s:%s", targetedGene.getGeneData().GeneName, targetedGeneRegion.getType().toString());

                if(targetedGeneRegion.useWholeRegion())
                {
                    panelRegion = new ProbeCandidate(region, ProbeSource.GENE, sourceInfo);
                }
                else
                {
                    ProbeCandidate probeCandidate = targetedGeneRegion.getSelectedProbe();

                    if(probeCandidate == null)
                    {
                        continue;
                    }

                    panelRegion = new ProbeCandidate(
                            region, ProbeSource.GENE, sourceInfo,
                            probeCandidate.getSequence(), probeCandidate.getGcContent(), probeCandidate.getQualityScore().get());
                }

                mPanelCache.addRegion(panelRegion);

                targetedGeneRegions.add(targetedGeneRegion);
            }
        }

        // write the outputs
        DataWriter.writeCandidates(mConfig.formOutputFilename(CANDIDATE_FILE_EXTENSION), targetedGeneRegions);
        DataWriter.writeTargetedGeneRegions(mConfig.formOutputFilename(GENE_REGION_FILE_EXTENSION), targetedGeneRegions);
    }
}
