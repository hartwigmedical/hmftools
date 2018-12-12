package com.hartwig.hmftools.svanalysis.annotators;

import java.util.List;
import java.util.Map;

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

    public DriverGeneAnnotator(DatabaseAccess dbAccess)
    {
        mDbAccess = dbAccess;
        mDriverCatalog = Lists.newArrayList();
        mSampleLOHData = null;

        initialiseGeneData();
    }

    public final List<DriverCatalog> getDriverCatalog() { return mDriverCatalog; }

    public void setSampleLohData(final Map<String, List<SvLOH>> data) { mSampleLohMap = data; }

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
        mSampleLOHData = mSampleLohMap.get(sampleId);

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
                annotateDeleteEvent(driverGene, region, false);
            else if(driverGene.driver() == DriverType.BIALLELIC)
                annotateDeleteEvent(driverGene, region, true);
            else if(driverGene.driver() == DriverType.AMP)
                annotateAmplification(driverGene);
        }
    }

    private void annotateDeleteEvent(final DriverCatalog driverGene, HmfTranscriptRegion region, boolean isSingleEvent)
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
            if (breakend.position() < region.start() && nextBreakend != null && nextBreakend.position() > region.start())
            {
                // take the lower of the copy numbers
                if (breakend.getCopyNumber(false) < nextBreakend.getCopyNumber(false))
                    minBreakend = breakend;
                else
                    minBreakend = nextBreakend;
            }
            else if (breakend.position() > region.end())
            {
                if(minBreakend != null)
                    break;

                // some event beyond the gene caused loss so need to continue on until it's found

                // for look for an LOH event starting from here
                for (int j = i; j < breakendList.size(); ++j)
                {
                    final SvBreakend lohBreakend = breakendList.get(j);
                    if (isLOHEvent(lohBreakend))
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
            LOGGER.warn("sample({}) driver gene({}) not allocated to SVs", mSampleId, driverGene.gene());
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

            setDriverGene(driverGene, startBreakend, "single-event start");
        }

        if(endBreakend != null)
            setDriverGene(driverGene, endBreakend, "single-event end");

        if(isSingleEvent)
            return;

        // look to next event
        final SvBreakend preStartBreakend = startBreakend != null ? findDeletionBreakend(breakendList, startBreakend.getChrPosIndex(), false, true) : null;
        final SvBreakend postEndBreakend = endBreakend != null ? findDeletionBreakend(breakendList, endBreakend .getChrPosIndex(), true, true) : null;

        if(preStartBreakend != null)
            setDriverGene(driverGene, preStartBreakend, "del-event pre-start");

        if(postEndBreakend != null)
            setDriverGene(driverGene, postEndBreakend, "del-event post-end");
    }

    private void setDriverGene(final DriverCatalog driverGene, final SvBreakend breakend, final String desc)
    {
        LOGGER.debug(String.format("gene(%s) %s breakend: %s cn(%.2f) cnChg(%.2f)",
                driverGene.gene(), desc, breakend.toString(), breakend.getCopyNumber(false),
                breakend.getSV().copyNumberChange(breakend.usesStart())));

        breakend.getSV().setDriveGene(driverGene, breakend.usesStart());
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
                    if(isLOHEvent(breakend))
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

    private boolean isLOHEvent(final SvBreakend breakend)
    {
        if (mSampleLOHData == null || mSampleLOHData.isEmpty())
            return false;

        for (final SvLOH lohEvent : mSampleLOHData)
        {
            if (lohEvent.StartSV.equals(breakend.getSV().id()) || lohEvent.EndSV.equals(breakend.getSV().id()))
                return true;
        }

        return false;
    }

    private void annotateAmplification(final DriverCatalog driverGene)
    {
        // find the cause - DUP, foldback, whole-chromatid duplication

    }

}
