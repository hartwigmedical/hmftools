package com.hartwig.hmftools.svanalysis.annotators;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_LOW_QUALITY;
import static com.hartwig.hmftools.svanalysis.types.SvLOH.LOH_NO_SV;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DriverGeneAnnotator
{
    private static final Logger LOGGER = LogManager.getLogger(DriverGeneAnnotator.class);

    final DatabaseAccess mDbAccess;

    private Map<String, HmfTranscriptRegion> mAllGenesMap;
    private List<DriverCatalog> mDriverCatalog;

    // references only
    private String mSampleId;
    private List<SvCluster> mClusters;
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private Map<String, List<SvLOH>> mSampleLohMap;
    private List<SvLOH> mSampleLOHData;
    private Map<String, Double> mChromosomeCopyNumberMap;

    public DriverGeneAnnotator(DatabaseAccess dbAccess)
    {
        mDbAccess = dbAccess;
        mDriverCatalog = Lists.newArrayList();
        mSampleLOHData = Lists.newArrayList();
        mChromosomeCopyNumberMap = null;

        initialiseGeneData();
    }

    public final List<DriverCatalog> getDriverCatalog() { return mDriverCatalog; }

    public void setChromosomeData(final Map<String, List<SvLOH>> lohData, Map<String, Double> chrCopyNumberMap)
    {
        mSampleLohMap = lohData;
        mChromosomeCopyNumberMap = chrCopyNumberMap;
    }

    private void initialiseGeneData()
    {
        SortedSetMultimap<String, HmfTranscriptRegion> genesByChromosomeMap = HmfGenePanelSupplier.allGenesPerChromosomeMap37();

        mAllGenesMap = Maps.newHashMap();
        for (final HmfTranscriptRegion region : genesByChromosomeMap.values())
        {
            mAllGenesMap.put(region.gene(), region);
        }
    }

    private final HmfTranscriptRegion findRegion(final String geneName)
    {
        return mAllGenesMap.get(geneName);
    }

    private void loadDriverCatalog(final String sampleId)
    {
        mDriverCatalog.clear();
        mDriverCatalog.addAll(mDbAccess.readDriverCatalog(sampleId));

        LOGGER.debug("sample({}) retrieved {} driver gene records", sampleId, mDriverCatalog.size());
    }

    public void annotateSVs(final String sampleId, final List<SvCluster> clusters, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        loadDriverCatalog(sampleId);

        if(mDriverCatalog.isEmpty())
            return;

        mSampleId = sampleId;
        mChrBreakendMap = chrBreakendMap;
        mClusters = clusters;

        mSampleLOHData.clear();
        List<SvLOH> sampleLohEvents = mSampleLohMap.get(sampleId);

        if(sampleLohEvents != null)
            mSampleLOHData.addAll(sampleLohEvents.stream().filter(x -> !x.Skipped).collect(Collectors.toList()));

        // Handle each of the 3 applicable types: DEL, BIALLELIC and AMP
        for(final DriverCatalog driverGene : mDriverCatalog)
        {
            HmfTranscriptRegion region = findRegion(driverGene.gene());

            if(region == null)
            {
                LOGGER.warn("driver gene({}) not found in all-genes map", driverGene.gene());
                continue;
            }

            if(driverGene.driver() == DriverType.DEL)
            {
                annotateDeleteEvent(driverGene, region);
            }
            else if(driverGene.driver() == DriverType.BIALLELIC)
            {
                annotateBiallelicEvent(driverGene, region);
            }
            else if(driverGene.driver() == DriverType.AMP)
            {
                annotateAmplification(driverGene, region);
            }
        }
    }

    private void annotateDeleteEvent(final DriverCatalog driverGene, HmfTranscriptRegion region)
    {
        /* DEL identification:
            - 1 or 2 SVs which caused this, start from DEL region (ie gene) and work out
            - gene copy number table - minCopyRegion < 0.5 makes it a DEL-type driver candidate
            - walk out from here to end of LOH event ie where (1- BAF*CN) > 0.5 again
            - see COLO829T and PTEN as an example
            - here there is a simple DEL and then most of the other chromatid has been lost
         */

        // first find the min copy number within the gene region
        // then walk out in both directions to find the SV which caused the loss
        // and then walk out again until heterozygosity is gained

        final List<SvBreakend> breakendList = mChrBreakendMap.get(region.chromosome());

        if (breakendList == null || breakendList.isEmpty())
            return;

        SvBreakend minBreakend = null; // breakend with lowest copy number covering any part of the gene region
        boolean isStartBreakend = true;

        for (int i = 0; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);
            final SvBreakend nextBreakend = i < breakendList.size() - 1 ? breakendList.get(i + 1) : null;

            // find the 2 breakends straddling the start of the gene region
            if ((breakend.position() < region.start() && nextBreakend != null && nextBreakend.position() > region.start())
            || (minBreakend == null && breakend.position() > region.start() && breakend.position() < region.end()))
            {
                // take the lower of the copy numbers
                if (nextBreakend == null || breakend.getCopyNumber(false) < nextBreakend.getCopyNumber(false))
                    minBreakend = breakend;
                else
                    minBreakend = nextBreakend;
            }
            else if (breakend.position() >= region.end())
            {
                if(minBreakend != null)
                    break;

                // some event beyond the gene caused loss so need to continue on until it's found

                // for look for an LOH event starting from here
                for (int j = i; j < breakendList.size(); ++j)
                {
                    final SvBreakend lohBreakend = breakendList.get(j);
                    if (isLOHEvent(lohBreakend, false))
                    {
                        minBreakend = lohBreakend;
                        isStartBreakend = false;
                        break;
                    }
                }

                if(minBreakend != null)
                    break;

                // otherwise select the first breakend facing away
                if(breakend.orientation() == -1)
                {
                    minBreakend = breakend;
                    break;
                }
            }
            else if (breakend.position() > region.start())
            {
                // another breakend inside the gene region
                if (breakend.getCopyNumber(false) < minBreakend.getCopyNumber(false))
                    minBreakend = breakend;
            }
        }

        if(minBreakend == null)
        {
            LOGGER.warn("sample({}) gene({}) not allocated to SVs", mSampleId, geneToStr(driverGene, region));
            return;
        }

        LOGGER.debug(String.format("gene(%s) min copy number at breakend(%s cn=%.2f cnChg=%.2f)",
                driverGene.gene(), minBreakend.toString(), minBreakend.getCopyNumber(false),
                minBreakend.getSV().copyNumberChange(minBreakend.usesStart())));

        SvBreakend startBreakend = isStartBreakend ? minBreakend : null;
        SvBreakend endBreakend = !isStartBreakend ? minBreakend : null;

        if(startBreakend != null)
        {
            // now walk forwards and backwards to find the next SVs
            int startIndex = startBreakend.getChrPosIndex();
            endBreakend = findDeletionBreakend(breakendList, startIndex, true, false);

            annotateDelSV(startBreakend, driverGene, region, "MIN");
        }

        if(endBreakend != null)
            annotateDelSV(endBreakend, driverGene, region, "MIN");

        // look to next event
        final SvBreakend preStartBreakend = startBreakend != null ? findDeletionBreakend(breakendList, startBreakend.getChrPosIndex(), false, true) : null;
        final SvBreakend postEndBreakend = endBreakend != null ? findDeletionBreakend(breakendList, endBreakend .getChrPosIndex(), true, true) : null;

        if(preStartBreakend != null)
            annotateDelSV(preStartBreakend, driverGene, region, "LOH");

        if(postEndBreakend != null)
            annotateDelSV(postEndBreakend, driverGene, region, "LOH");
    }

    private void annotateBiallelicEvent(final DriverCatalog driverGene, HmfTranscriptRegion region)
    {
        // for biallelic events, find the straddling LOH event
        if(mSampleLOHData.isEmpty())
            return;

        final List<SvBreakend> breakendList = mChrBreakendMap.get(region.chromosome());

        if (breakendList == null || breakendList.isEmpty())
            return;

        // find any LOH which cross over all or a part of this gene region
        for (final SvLOH lohEvent : mSampleLOHData)
        {
            if(lohEvent.PosStart > region.end() || lohEvent.PosEnd < region.start())
                continue;

            // now find the corresponding breakends
            SvBreakend startBreakend = null;
            SvBreakend endBreakend = null;

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);

                if (lohEvent.StartSV.equals(breakend.getSV().id()))
                    startBreakend = breakend;

                if (lohEvent.EndSV.equals(breakend.getSV().id()))
                    endBreakend = breakend;

                if ((startBreakend != null || lohEvent.StartSV.equals(LOH_NO_SV))
                && (endBreakend != null || lohEvent.EndSV.equals(LOH_NO_SV)))
                {
                    break;
                }
            }

            if(startBreakend == null && endBreakend == null)
            {
                LOGGER.warn("sample({}) gene({}) not allocated to SVs",
                        mSampleId, geneToStr(driverGene, region));
                return;
            }

            if(startBreakend != null)
                annotateDelSV(startBreakend, driverGene, region, "LOH");

            if(endBreakend != null)
                annotateDelSV(endBreakend, driverGene, region, "LOH");

            break;
        }
    }

    private SvBreakend findDeletionBreakend(final List<SvBreakend> breakendList, int startIndex, boolean walkForwards, boolean requireGOH)
    {
        int index = walkForwards ? startIndex + 1 : startIndex - 1;

        while(index >= 0 && index <= breakendList.size() - 1)
        {
            final SvBreakend breakend = breakendList.get(index);

            // this first variant following the DEL region ought to be facing away
            if ((walkForwards && breakend.orientation() == -1) || (!walkForwards && breakend.orientation() == 1))
            {
                if(requireGOH)
                {
                    if(isLOHEvent(breakend, !walkForwards))
                        return breakend;
                }
                else
                {
                    return breakend;
                }
            }

            index = walkForwards ? index + 1 : index - 1;
        }

        return null;
    }

    private boolean isLOHEvent(final SvBreakend breakend, boolean checkStart)
    {
        for (final SvLOH lohEvent : mSampleLOHData)
        {
            if ((checkStart && lohEvent.StartSV.equals(breakend.getSV().id()))
            || (!checkStart && lohEvent.EndSV.equals(breakend.getSV().id())))
            {
                return true;
            }
        }

        return false;
    }

    private void annotateAmplification(final DriverCatalog driverGene, HmfTranscriptRegion region)
    {
        // find the cause - DUP, foldback, otherwise assume whole-chromatid duplication

        final List<SvBreakend> breakendList = mChrBreakendMap.get(region.chromosome());

        if (breakendList == null || breakendList.isEmpty())
            return;

        // trying to find the breakends which amplify this gene
        // take any INV, DUP or TI which straddles the gene
        int candidateCount = 0;

        for (int i = 0; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);

            if(breakend.position() > region.start())
                break;

            final SvVarData varStart = breakend.getSV();

            if(varStart.getCluster().getResolvedType() == RESOLVED_TYPE_LOW_QUALITY)
                continue;

            if(varStart.type() == DUP || varStart.type() == INV)
            {
                if(varStart.position(true) <= region.start() && varStart.position(false) >= region.end())
                {
                    LOGGER.info(String.format("sample(%s) cluster(%s) gene(%s) single SV(%s %s) cn(%.2f) cnChg(%.2f)",
                            mSampleId, varStart.getCluster().id(), geneToStr(driverGene, region), varStart.posId(), varStart.type(),
                            varStart.copyNumber(true), varStart.copyNumberChange(true)));

                    annotateSV(varStart, driverGene, true, "SV");
                    annotateSV(varStart, driverGene, false, "SV");
                    ++candidateCount;
                    continue;
                }
            }

            final SvLinkedPair tiPair = varStart.getLinkedPair(breakend.usesStart());

            if(tiPair == null)
                continue;

            if(breakend.position() + tiPair.length() < region.end())
                continue;

            final SvVarData varEnd = tiPair.getOtherSV(varStart);
            boolean v1Start = varStart == tiPair.first() ? tiPair.firstLinkOnStart() : tiPair.secondLinkOnStart();
            boolean v2Start = varEnd == tiPair.first() ? tiPair.firstLinkOnStart() : tiPair.secondLinkOnStart();

            LOGGER.info(String.format("sample(%s) cluster(%d fb=%s) gene(%s) SVs start(%s cn=%.2f cnChg=%.2f) end(%s cn=%.2f cnChg=%.2f) in linked pair",
                    mSampleId, varStart.getCluster().id(), varStart.getCluster().getFoldbacks().size(), geneToStr(driverGene, region),
                    varStart.posId(), varStart.copyNumber(v1Start), varStart.copyNumberChange(v1Start),
                    varEnd.posId(), varEnd.copyNumber(v2Start), varEnd.copyNumberChange(v2Start)));

            annotateSV(varStart, driverGene, v1Start, "TI");
            annotateSV(varEnd, driverGene, v2Start, "TI");
            ++candidateCount;
        }

        if(candidateCount == 0)
        {
            // other likely cause is whole-chromatid amplification
            Double chrCopyNumber = mChromosomeCopyNumberMap.get(region.chromosome());

            if (chrCopyNumber != null)
            {
                LOGGER.info(String.format("sample(%s) gene(AMP: %s) chr(%s) copy-number(%.2f)",
                        mSampleId, driverGene.gene(), region.chromosome(), chrCopyNumber));
            }
        }
    }


    private void annotateDelSV(final SvBreakend breakend, DriverCatalog driverGene, HmfTranscriptRegion region, final String desc)
    {
        final SvVarData var = breakend.getSV();

        LOGGER.info(String.format("sample(%s) cluster(%d) gene(%s) single SV(%s %s) cn(%.2f) cnChg(%.2f) as %s",
                mSampleId, var.getCluster().id(), geneToStr(driverGene, region), var.posId(), var.type(),
                var.copyNumber(breakend.usesStart()), var.copyNumberChange(breakend.usesStart()), desc));

        annotateSV(var, driverGene, breakend.usesStart(), desc);
    }

    private static final String geneToStr(DriverCatalog driverGene, HmfTranscriptRegion region)
    {
        return String.format("%s: %s %s:%d-%d",
                driverGene.driver(), driverGene.gene(), region.chromosome(), region.start(), region.end());
    }

    private static void annotateSV(final SvVarData var, DriverCatalog driverGene, boolean isStart, final String desc)
    {
        var.setDriveGene(String.format("%s;%s;%s", driverGene.driver(), driverGene.gene(), desc), isStart);
    }

}
