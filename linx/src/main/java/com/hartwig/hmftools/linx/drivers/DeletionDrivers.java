package com.hartwig.hmftools.linx.drivers;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.MAX_COPY_NUM_DIFF;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.MAX_COPY_NUM_DIFF_PERC;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.TOTAL_CN_LOSS;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.expectSingleChromosome;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.HOM_DEL_DISRUPTION;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.HOM_DUP_DISRUPTION;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_ARM;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_CHR;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_SV_CENTRO;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_SV_TELO;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_ARM_SV;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_DEL;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.SvBreakend;

public class DeletionDrivers
{
    private final DriverDataCache mDataCache;
    private final Map<String,List<String>> mReportableDisruptionGeneTranscripts;
    private final boolean mHomDisAllGenes;

    public DeletionDrivers(final Map<String,List<String>> disruptionGeneTranscripts, final DriverDataCache dataCache, boolean homDisAllGenes)
    {
        mDataCache = dataCache;
        mReportableDisruptionGeneTranscripts = disruptionGeneTranscripts;
        mHomDisAllGenes = homDisAllGenes;
    }

    public void annotateDeleteEvent(final DriverGeneData dgData)
    {
        // skip if a germline event causes this DEL or the calculated min CN doesn't support it
        if(dgData.CopyNumberRegion.IsGermline || dgData.CopyNumberRegion.MinCopyNumber >= 0.51)
            return;

        // find the LOH and hom-loss events which caused this DEL
        int minRegionStart = dgData.CopyNumberRegion.RegionStart;
        int minRegionEnd = dgData.CopyNumberRegion.RegionEnd;

        if(expectSingleChromosome(mDataCache.isMale(), dgData.GeneData.Chromosome))
        {
            for(final HomLossEvent homLoss : mDataCache.CopyNumberData.getHomLossData())
            {
                // allow there to be more than one?
                if(homLoss.PosStart <= minRegionEnd && homLoss.PosEnd >= minRegionStart)
                {
                    LNX_LOGGER.debug("gene({}) minCnRegion({} -> {}) covered by hom-loss({})",
                            dgData.GeneData.GeneName, minRegionStart, minRegionEnd, homLoss);

                    // DEL covers the whole gene or extends to end of arm
                    DriverGeneEvent event = new DriverGeneEvent(DriverEventType.DEL);
                    event.setHomLossEvent(homLoss);

                    // one or both breakends can be null if not matched
                    event.addSvBreakendPair(homLoss.getBreakend(true), homLoss.getBreakend(false), SV_DRIVER_TYPE_DEL);
                    dgData.addEvent(event);
                    break;
                }
            }

            return;
        }

        for(final LohEvent lohEvent : mDataCache.CopyNumberData.getLohData())
        {
            if(!lohEvent.Chromosome.equals(dgData.GeneData.Chromosome))
                continue;

            // the LOH needs to straddle one or the other of the min-gene breakends
            if(lohEvent.PosStart > minRegionEnd || lohEvent.PosEnd < minRegionStart)
                continue;

            LNX_LOGGER.debug("gene({}) minCnRegion({} -> {}) covered by LOH({})",
                    dgData.GeneData.GeneName, minRegionStart, minRegionEnd, lohEvent);

            if((lohEvent.PosStart < minRegionEnd && lohEvent.PosEnd > minRegionStart) || !lohEvent.doubleSvEvent())
            {
                // LOH covers the whole gene or extends to end of arm (ie not just 2 SVs
                if(lohEvent.doubleSvEvent())
                {
                    DriverGeneEvent event = new DriverGeneEvent(LOH);
                    event.setLohEvent(lohEvent);

                    // one or both breakends can be null if not matched
                    event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_DEL);
                    dgData.addEvent(event);
                }
                else if(lohEvent.armLoss())
                {
                    dgData.addEvent(new DriverGeneEvent(LOH_ARM));
                }
                else if(lohEvent.chromosomeLoss())
                {
                    dgData.addEvent(new DriverGeneEvent(LOH_CHR));
                }
                else if(lohEvent.isSvEvent())
                {
                    // call SV + rest of arm loss an LOH as well
                    DriverEventType eventType = lohEvent.telomereLoss() ? LOH_SV_TELO : LOH_SV_CENTRO;
                    DriverGeneEvent event = new DriverGeneEvent(eventType);
                    event.setLohEvent(lohEvent);
                    event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_ARM_SV);
                    dgData.addEvent(event);
                }
            }
            else
            {
                // cannot assign one event as the LOH and one as the DEL, so assign both as DELs
                DriverGeneEvent event = new DriverGeneEvent(DriverEventType.DEL);
                event.setLohEvent(lohEvent);

                // one or both breakends can be null if not matched
                event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_DEL);
                dgData.addEvent(event);
            }

            for(final HomLossEvent homLoss : lohEvent.getHomLossEvents())
            {
                // allow there to be more than one?
                if(homLoss.PosStart <= minRegionEnd && homLoss.PosEnd >= minRegionStart)
                {
                    LNX_LOGGER.debug("gene({}) minCnRegion({} -> {}) covered by hom-loss({})",
                            dgData.GeneData.GeneName, minRegionStart, minRegionEnd, homLoss);

                    // DEL covers the whole gene or extends to end of arm
                    DriverGeneEvent event = new DriverGeneEvent(DriverEventType.DEL);
                    event.setHomLossEvent(homLoss);

                    // one or both breakends can be null if not matched
                    event.addSvBreakendPair(homLoss.getBreakend(true), homLoss.getBreakend(false), SV_DRIVER_TYPE_DEL);
                    dgData.addEvent(event);
                    break;
                }
            }

            break;
        }
    }

    public void annotateBiallelicEvent(final DriverGeneData dgData)
    {
        // look for an LOH covering any part of the coding region
        if(dgData.TransData.CodingStart == null || dgData.TransData.CodingEnd == null)
            return;

        int codingStart = dgData.TransData.CodingStart;
        int codingEnd = dgData.TransData.CodingEnd;

        for(final LohEvent lohEvent : mDataCache.CopyNumberData.getLohData())
        {
            if(!lohEvent.Chromosome.equals(dgData.GeneData.Chromosome))
                continue;

            if(lohEvent.PosStart > codingEnd || lohEvent.PosEnd < codingStart)
                continue;

            LNX_LOGGER.debug("gene({}) coding region({} -> {}) covered by LOH({})",
                    dgData.GeneData.GeneName, codingStart, codingEnd, lohEvent);

            if(lohEvent.doubleSvEvent())
            {
                DriverGeneEvent event = new DriverGeneEvent(LOH);
                event.setLohEvent(lohEvent);

                // one or both breakends can be null if not matched
                event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_DEL);
                dgData.addEvent(event);
            }
            else if(lohEvent.armLoss())
            {
                dgData.addEvent(new DriverGeneEvent(LOH_ARM));
            }
            else if(lohEvent.chromosomeLoss())
            {
                dgData.addEvent(new DriverGeneEvent(LOH_CHR));
            }
            else if(lohEvent.isSvEvent())
            {
                DriverEventType eventType = lohEvent.telomereLoss() ? LOH_SV_TELO : LOH_SV_CENTRO;
                DriverGeneEvent event = new DriverGeneEvent(eventType);
                event.setLohEvent(lohEvent);
                event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_ARM_SV);
                dgData.addEvent(event);
            }

            break;
        }
    }

    public List<DriverGeneData> findDisruptiveDelDrivers(final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        // find homozygous disruptions by looking for either:
        // 1. Evidence of copy number dropping to zero within a gene from disruptive SVs other than simple (unphased) DELs
        // 2. Duplication of exons from DUPs or other complex clusters
        final List<String> delDriverGeneIds = mDataCache.getDriverGeneDataList().stream()
                .filter(x -> x.DriverData.driver() == DriverType.DEL)
                .map(x -> x.GeneData.GeneId).collect(Collectors.toList());

        final List<DriverGeneData> disDelDrivers = Lists.newArrayList();

        for(Map.Entry<String,List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            for(final SvBreakend breakend : entry.getValue())
            {
                boolean delType = breakend.orientation() == 1 && breakend.getDBLink() != null;
                boolean dupType = breakend.orientation() == -1 && breakend.type() == DUP;

                if(!delType && !dupType)
                    continue;

                final List<BreakendGeneData> genesList = breakend.getSV().getGenesList(breakend.usesStart()).stream()
                        .filter(x -> !delDriverGeneIds.contains(x.StableId))
                        .filter(x -> mHomDisAllGenes || mReportableDisruptionGeneTranscripts.containsKey(x.StableId))
                        .collect(Collectors.toList());

                if(genesList.isEmpty())
                    continue;

                // consider any disruptive canonical or other configured reportable transcript
                for(final BreakendGeneData gene : genesList)
                {
                    for(BreakendTransData transcript : gene.transcripts())
                    {
                        if(transcript == null || !transcript.isDisruptive())
                            continue;

                        if(!transcript.isCanonical())
                        {
                            if(!mReportableDisruptionGeneTranscripts.get(gene.StableId).contains(transcript.transName()))
                                continue;
                        }

                        if(delType)
                        {
                            checkDelDisruption(breakend, transcript, delDriverGeneIds, disDelDrivers);
                        }
                        else
                        {
                            checkDupDisruption(breakend, transcript, delDriverGeneIds, disDelDrivers);
                        }
                    }
                }
            }
        }

        return disDelDrivers;
    }

    private void checkDelDisruption(
            final SvBreakend breakend, final BreakendTransData transcript,
            final List<String> delDriverGeneIds, final List<DriverGeneData> disDelDrivers)
    {
        final DbPair dbLink = breakend.getDBLink();

        // calculate the copy number over the deletion bridge section
        double cnLowSide = breakend.copyNumberLowSide();

        double otherSvJcn = dbLink.getOtherBreakend(breakend).jcn();

        // account for an overlapping DB by subtracting the ploidy of the overlap
        if(dbLink.length() < 0)
            cnLowSide -= otherSvJcn;

        if(cnLowSide >= TOTAL_CN_LOSS)
            return;

        // both the SVs involved in the deletion must be disruptive, ie cannot be simple intronic DELs
        final String geneId = transcript.gene().StableId;

        final SvBreakend otherBreakend = dbLink.getOtherBreakend(breakend);

        final BreakendGeneData matchingGene = otherBreakend.getSV().getGenesList(otherBreakend.usesStart()).stream()
                .filter(x -> x.StableId.equals(geneId))
                .findFirst().orElse(null);

        if(matchingGene == null)
            return;

        // the other transcript must also be disruptive
        if(matchingGene.transcripts().stream().noneMatch(x -> x.transName().equals(transcript.transName()) && x.isDisruptive()))
            return;

        delDriverGeneIds.add(geneId);

        LNX_LOGGER.debug("breakend({}) cluster({}) geneTrans({}:{}) cause homozygous disruption for cnLowSide({}) dbLength({}) otherSvJcn({})",
                breakend, breakend.getCluster().id(), transcript.geneName(), transcript.transName(),
                formatJcn(cnLowSide), dbLink.length(), formatJcn(otherSvJcn));

        DriverGeneData dgData = mDataCache.createDriverData(transcript.gene(), transcript.TransData);

        if(dgData != null)
        {
            DriverGeneEvent event = new DriverGeneEvent(HOM_DEL_DISRUPTION);
            event.addSvBreakendPair(breakend, otherBreakend, "DB");
            event.setCluster(breakend.getCluster());
            dgData.addEvent(event);
            disDelDrivers.add(dgData);
        }
    }

    private void checkDupDisruption(
            final SvBreakend breakend, final BreakendTransData transcript,
            final List<String> delDriverGeneIds, final List<DriverGeneData> disDelDrivers)
    {
        final String geneId = transcript.gene().StableId;

        // DUP must be wholly contained within the same gene
        if(breakend.getSV().getGenesList(!breakend.usesStart()).stream().noneMatch(x -> x.StableId.equals(geneId)))
            return;

        final SvBreakend otherBreakend = breakend.getOtherBreakend();

        final BreakendGeneData matchingGene = otherBreakend.getSV().getGenesList(otherBreakend.usesStart()).stream()
                .filter(x -> x.StableId.equals(geneId))
                .findFirst().orElse(null);

        if(matchingGene == null)
            return;

        // the other transcript must also be disruptive
        if(matchingGene.transcripts().stream().noneMatch(x -> x.transName().equals(transcript.transName()) && x.isDisruptive()))
            return;

        double cnLowSideStart = breakend.copyNumberLowSide();
        double cnLowSideEnd = otherBreakend.copyNumberLowSide();
        double jcn = breakend.jcn();
        double jcnThreshold = max(jcn * (1 + MAX_COPY_NUM_DIFF_PERC), jcn + MAX_COPY_NUM_DIFF);

        if(cnLowSideStart < jcnThreshold && cnLowSideEnd < jcnThreshold)
        {
            delDriverGeneIds.add(geneId);

            LNX_LOGGER.debug("DUP({}) cluster({}) gene({}:{}) cause homozygous disruption cnLowSide({} & {}) jcn({})",
                    breakend, breakend.getCluster().id(), transcript.geneName(), transcript.transName(),
                    transcript.geneName(), breakend.getCluster().id(), breakend.getSV().id(),
                    formatJcn(cnLowSideStart), formatJcn(cnLowSideEnd), formatJcn(jcn));

            DriverGeneData dgData = mDataCache.createDriverData(transcript.gene(), transcript.TransData);

            if(dgData != null)
            {
                DriverGeneEvent event = new DriverGeneEvent(HOM_DUP_DISRUPTION);
                event.addSvBreakendPair(breakend, otherBreakend, "DUP");
                event.setCluster(breakend.getCluster());
                dgData.addEvent(event);
                disDelDrivers.add(dgData);
            }
        }
    }
}
