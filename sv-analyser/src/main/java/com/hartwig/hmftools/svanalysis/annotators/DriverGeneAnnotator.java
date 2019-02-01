package com.hartwig.hmftools.svanalysis.annotators;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.CENTROMERE_CN;
import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.P_ARM_TELOMERE_CN;
import static com.hartwig.hmftools.svanalysis.analysis.CNAnalyser.Q_ARM_TELOMERE_CN;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser.SHORT_TI_LENGTH;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_LOW_QUALITY;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvLOH.LOH_NO_SV;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvLOH;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DriverGeneAnnotator
{
    private static final Logger LOGGER = LogManager.getLogger(DriverGeneAnnotator.class);

    final DatabaseAccess mDbAccess;
    SvGeneTranscriptCollection mGeneTranscriptCollection;

    private Map<String, HmfTranscriptRegion> mAllGenesMap;
    private List<DriverCatalog> mDriverCatalog;
    private List<GeneCopyNumber> mGeneCopyNumberData;
    private BufferedWriter mFileWriter;
    private String mOutputDir;
    private double mSamplePloidy;

    // references only
    private String mSampleId;
    private List<SvCluster> mClusters;
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private Map<String, List<SvLOH>> mSampleLohMap;
    private List<SvLOH> mSampleLOHData;
    private Map<String, double[]> mChrCopyNumberMap;

    public DriverGeneAnnotator(DatabaseAccess dbAccess, SvGeneTranscriptCollection geneTranscriptCollection, final String outputDir)
    {
        mDbAccess = dbAccess;
        mGeneTranscriptCollection = geneTranscriptCollection;

        mDriverCatalog = Lists.newArrayList();
        mGeneCopyNumberData = Lists.newArrayList();
        mSampleLOHData = Lists.newArrayList();
        mChrCopyNumberMap = null;
        mFileWriter = null;
        mOutputDir = outputDir;
        mSamplePloidy = 0;

        initialiseGeneData();
    }

    public void setLohData(final Map<String, List<SvLOH>> lohData) { mSampleLohMap = lohData; }
    public void setChrCopyNumberMap(Map<String, double[]> chrCopyNumberMap) { mChrCopyNumberMap = chrCopyNumberMap; }
    public void setSamplePloidy(double ploidy)
    {
        mSamplePloidy = ploidy;
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

    private void loadGeneCopyNumberData(final String sampleId)
    {
        mGeneCopyNumberData.clear();
        mGeneCopyNumberData = mDbAccess.readGeneCopynumbers(sampleId);
    }

    public void annotateSVs(final String sampleId, final List<SvCluster> clusters, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mSampleId = sampleId;
        mChrBreakendMap = chrBreakendMap;
        mClusters = clusters;

        checkPseudoGeneAnnotations();

        loadDriverCatalog(sampleId);

        if (mDriverCatalog.isEmpty())
            return;

        // loadGeneCopyNumberData(sampleId);

        mSampleLOHData.clear();
        List<SvLOH> sampleLohEvents = mSampleLohMap.get(sampleId);

        if (sampleLohEvents != null)
            mSampleLOHData.addAll(sampleLohEvents.stream().filter(x -> !x.Skipped).collect(Collectors.toList()));

        // Handle each of the 3 applicable types: DEL, BIALLELIC and AMP
        for (final DriverCatalog driverGene : mDriverCatalog)
        {
            HmfTranscriptRegion region = findRegion(driverGene.gene());

            if (region == null)
            {
                LOGGER.warn("driver gene({}) not found in all-genes map", driverGene.gene());
                continue;
            }

            final List<SvBreakend> breakendList = mChrBreakendMap.get(region.chromosome());

            if (breakendList == null || breakendList.isEmpty())
            {
                writeDriverData(driverGene, region);
                continue;
            }

            if (driverGene.driver() == DriverType.DEL)
            {
                annotateDeleteEvent(driverGene, region, breakendList);
            }
            else if (driverGene.driver() == DriverType.AMP)
            {
                annotateAmplification(driverGene, region, breakendList);
            }
            else
            {
                // treat DNDS and HOTSPOT as potentially biallelic events
                annotateBiallelicEvent(driverGene, region, breakendList);
            }
        }
    }

    private void annotateDeleteEvent(final DriverCatalog driverGene, HmfTranscriptRegion region, final List<SvBreakend> breakendList)
    {
        /* DEL identification:
            - 1 or 2 SVs which caused this, start from DEL region (ie gene) and work out
            - gene copy number table - minCopyRegion < 0.5 makes it a DEL-type driver candidate
            - walk out from here to end of LOH event ie where (1- BAF*CN) > 0.5 again
            - see COLO829T and PTEN as an example - here there is a simple DEL and then most of the other chromatid has been lost
         */

        // first find the min copy number within the gene region
        // then walk out in both directions to find the SV which caused the loss
        // and then walk out again until heterozygosity is gained (ie the end of an LOH)

        SvBreakend minBreakend = null; // breakend with lowest copy number covering any part of the gene region
        boolean isStartBreakend = true;

        for (int i = 0; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);
            final SvBreakend nextBreakend = i < breakendList.size() - 1 ? breakendList.get(i + 1) : null;

            // find the 2 breakends straddling the start of the gene region, and compare their CNs in the -1 / upwards direction
            if ((breakend.position() < region.start() && nextBreakend != null && nextBreakend.position() > region.start())
            || (minBreakend == null && breakend.position() > region.start() && breakend.position() < region.end()))
            {
                // take the lower of the copy numbers (upwards-facing)
                if (nextBreakend == null || breakend.getCopyNumber(false) < nextBreakend.getCopyNumber(false)
                || nextBreakend.position() >= region.end())
                {
                    minBreakend = breakend;
                }
                else
                {
                    minBreakend = nextBreakend;
                }
            }
            else if(minBreakend == null && nextBreakend == null)
            {
                // no more breakends to consider - does this breakend definitely contribute to the DEL?
                minBreakend = breakend;
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
                // check for a lower CN from another breakend inside the gene region
                if (breakend.getCopyNumber(false) < minBreakend.getCopyNumber(false))
                    minBreakend = breakend;
            }
        }

        if(minBreakend == null)
        {
            LOGGER.debug("sample({}) gene({}) not allocated to SVs", mSampleId, geneToStr(driverGene, region));
            writeDriverData(driverGene, region);
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

        // find straddling LOH for this DEL
        SvBreakend preStartBreakend = null;
        SvBreakend postEndBreakend = null;

        SvLOH matchedLohEvent = null;
        if(startBreakend != null && endBreakend != null)
        {
            for (final SvLOH lohEvent : mSampleLOHData)
            {
                if(lohEvent.Chromosome.equals(startBreakend.chromosome())
                && lohEvent.PosStart <= startBreakend.position() && lohEvent.PosEnd >= endBreakend.position())
                {
                    matchedLohEvent = lohEvent;
                    preStartBreakend = lohEvent.getBreakend(true);
                    postEndBreakend = lohEvent.getBreakend(false);
                    break;
                }
            }
        }

        // final SvBreakend preStartBreakend = startBreakend != null ? findDeletionBreakend(breakendList, startBreakend.getChrPosIndex(), false, true) : null;
        // final SvBreakend postEndBreakend = endBreakend != null ? findDeletionBreakend(breakendList, endBreakend.getChrPosIndex(), true, true) : null;

        // look to next event
        if(preStartBreakend == null && postEndBreakend == null)
        {
            // search for a whole arm or chromatid LOH event
            writeDriverData(driverGene, region);
        }
        else if(preStartBreakend != null && postEndBreakend != null)
        {
            annotateDelSV(preStartBreakend, driverGene, region, "LOH");
            annotateDelSV(postEndBreakend, driverGene, region, "LOH");
        }
        else if(preStartBreakend != null)
        {
            annotateDelSV(preStartBreakend, driverGene, region, "LOH_" + matchedLohEvent.SegEnd);
        }
        else
        {
           annotateDelSV(postEndBreakend, driverGene, region, "LOH_" + matchedLohEvent.SegStart);
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

    private void annotateBiallelicEvent(final DriverCatalog driverGene, HmfTranscriptRegion region, final List<SvBreakend> breakendList)
    {
        // for biallelic events, find the straddling LOH event
        for (final SvLOH lohEvent : mSampleLOHData)
        {
            if(!lohEvent.Chromosome.equals(region.chromosome()))
                continue;

            if(lohEvent.PosStart > region.end() || lohEvent.PosEnd < region.start())
                continue;

            // now find the corresponding breakends
            SvBreakend startBreakend = lohEvent.getBreakend(true);
            SvBreakend endBreakend = lohEvent.getBreakend(false);

            if(startBreakend == null && endBreakend == null)
            {
                LOGGER.debug("sample({}) gene({}) LOH not found for driver gene",
                        mSampleId, geneToStr(driverGene, region));
                break;
            }


            if(startBreakend != null && endBreakend != null)
            {
                annotateDelSV(startBreakend, driverGene, region, "LOH");
                annotateDelSV(endBreakend, driverGene, region, "LOH");
            }
            else if(startBreakend != null)
            {
                annotateDelSV(startBreakend, driverGene, region, "LOH_" + lohEvent.SegEnd);
            }
            else
            {
                annotateDelSV(endBreakend, driverGene, region, "LOH_" + lohEvent.SegStart);
            }

            return;
        }

        writeDriverData(driverGene, region);
    }

    private void annotateDelSV(final SvBreakend breakend, DriverCatalog driverGene, HmfTranscriptRegion region, final String desc)
    {
        final SvVarData var = breakend.getSV();

        LOGGER.debug(String.format("sample(%s) cluster(%d) gene(%s) single SV(%s %s) cn(%.2f) cnChg(%.2f) as %s",
                mSampleId, var.getCluster().id(), geneToStr(driverGene, region), var.posId(), var.type(),
                var.copyNumber(breakend.usesStart()), var.copyNumberChange(breakend.usesStart()), desc));

        var.setDriveGene(String.format("%s;%s;%s", driverGene.driver(), driverGene.gene(), desc), breakend.usesStart());
        writeDriverData(driverGene, region, var, breakend.usesStart(), desc);
    }

    private void annotateAmplification(final DriverCatalog driverGene, HmfTranscriptRegion region, final List<SvBreakend> breakendList)
    {
        // find the cause - DUP, foldback, otherwise assume whole-chromatid duplication

        // trying to find the breakends which amplify this gene
        // take any INV, DUP or TI which straddles the gene
        double maxCopyNumber = 0;
        SvVarData maxSvStart = null;
        SvVarData maxSvEnd = null;
        SvBreakend foldbackBreakend = null;

        for (int i = 0; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);

            if(breakend.position() > region.start())
            {
                // past the gene but AMP may be explained by an unchained foldback
                if(foldbackBreakend == null && breakend.orientation() == 1
                && !breakend.getSV().getFoldbackLink(breakend.usesStart()).isEmpty())
                {
                    foldbackBreakend = breakend;
                    break;
                }

                continue; // looking for a foldback as an explanation
            }

            final SvVarData varStart = breakend.getSV();

            if(varStart.getCluster().getResolvedType() == RESOLVED_TYPE_LOW_QUALITY)
                continue;

            if(breakend.orientation() == -1 && !varStart.getFoldbackLink(breakend.usesStart()).isEmpty())
                foldbackBreakend = breakend;

            if(varStart.type() == DUP || varStart.type() == INV)
            {
                if(varStart.position(true) <= region.start() && varStart.position(false) >= region.end())
                {
                    LOGGER.debug(String.format("sample(%s) cluster(%s) gene(%s) single SV(%s %s) cn(%.2f) cnChg(%.2f)",
                            mSampleId, varStart.getCluster().id(), geneToStr(driverGene, region), varStart.posId(), varStart.type(),
                            varStart.copyNumber(true), varStart.copyNumberChange(true)));

                    annotateSV(varStart, driverGene,  region,  true, "SV");
                    annotateSV(varStart, driverGene, region, false, "SV");

                    if(varStart.copyNumberChange(true) > maxCopyNumber)
                    {
                        maxCopyNumber = varStart.copyNumberChange(true);
                        maxSvStart = varStart;
                        maxSvEnd = varStart;
                    }

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

            LOGGER.debug(String.format("sample(%s) cluster(%d fb=%s) gene(%s) SVs start(%s cn=%.2f cnChg=%.2f) end(%s cn=%.2f cnChg=%.2f) in linked pair",
                    mSampleId, varStart.getCluster().id(), varStart.getCluster().getFoldbacks().size(), geneToStr(driverGene, region),
                    varStart.posId(), varStart.copyNumber(v1Start), varStart.copyNumberChange(v1Start),
                    varEnd.posId(), varEnd.copyNumber(v2Start), varEnd.copyNumberChange(v2Start)));

            annotateSV(varStart, driverGene, region, v1Start, "TI");
            annotateSV(varEnd, driverGene, region, v2Start, "TI");

            if(varStart.copyNumberChange(true) > maxCopyNumber)
            {
                maxCopyNumber = varStart.copyNumberChange(true);
                maxSvStart = varStart;
                maxSvEnd = varEnd;
            }
        }

        if(maxSvStart != null && maxSvEnd != null)
            return;

        if(foldbackBreakend != null)
            annotateSV(foldbackBreakend.getSV(), driverGene, region, foldbackBreakend.usesStart(), "FB");
        else
            writeDriverData(driverGene, region);
    }

    private void annotateSV(final SvVarData var, final DriverCatalog driverGene, final HmfTranscriptRegion region, boolean isStart, final String desc)
    {
        var.setDriveGene(String.format("%s;%s;%s", driverGene.driver(), driverGene.gene(), desc), isStart);
        writeDriverData(driverGene, region, var, isStart, desc);
    }

    private static final String geneToStr(DriverCatalog driverGene, HmfTranscriptRegion region)
    {
        return String.format("%s: %s %s:%d-%d",
                driverGene.driver(), driverGene.gene(), region.chromosome(), region.start(), region.end());
    }

    private void writeDriverData(final DriverCatalog driverGene, final HmfTranscriptRegion region)
    {
        writeDriverData(driverGene, region, null, true, "");
    }

    private final SvLOH findLOHEventForRegion(final HmfTranscriptRegion region)
    {
        // check for any LOH not linked to an SV which straddles this region
        for (final SvLOH lohEvent : mSampleLOHData)
        {
            if (lohEvent.Chromosome.equals(region.chromosome())
            && lohEvent.PosStart <= region.start() && lohEvent.PosEnd >= region.end())
            {
                return lohEvent;
            }
        }

        return null;
    }

    private final GeneCopyNumber findGeneCopyNumber(final DriverCatalog driverGene)
    {
        for(GeneCopyNumber geneCN : mGeneCopyNumberData)
        {
            if(driverGene.gene().equals(geneCN.gene()))
                return geneCN;
        }

        return null;
    }

    private void writeDriverData(final DriverCatalog driverGene, final HmfTranscriptRegion region, final SvVarData var, boolean isStart,
            final String matchInfo)
    {
        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                if(!outputFileName.endsWith("/"))
                    outputFileName += File.separator;

                outputFileName += "SVA_DRIVERS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Gene,GeneType,DriverType,ClusterId,SvId,IsStart,MatchInfo");
                mFileWriter.write(",SamplePloidy,Chromosome,Arm,MinCN,CentromereCN,TelomereCN");
                mFileWriter.newLine();
            }

            BufferedWriter writer = mFileWriter;

            writer.write(String.format("%s,%s,%s,%s", mSampleId, driverGene.gene(), driverGene.category(), driverGene.driver()));

            if(var != null)
            {
                writer.write(String.format(",%d,%s,%s,%s", var.getCluster().id(), var.origId(), isStart, matchInfo));
            }
            else
            {
                final SvLOH lohEvent = driverGene.driver() != DriverType.AMP ? findLOHEventForRegion(region) : null;

                if (lohEvent != null)
                        writer.write(String.format(",-1,,,%s;%s", lohEvent.SegStart, lohEvent.SegEnd));
                    else
                        writer.write(String.format(",-1,,,"));
            }

            double[] cnData = mChrCopyNumberMap.get(region.chromosome());

            GeneCopyNumber geneCN = findGeneCopyNumber(driverGene);

            final String arm = var != null ? var.arm(isStart) : getChromosomalArm(region.chromosome(), region.geneStart());

            writer.write(String.format(",%.2f,%s,%s,%.2f,%.2f,%.2f",
                    mSamplePloidy, region.chromosome(), arm, geneCN != null ? geneCN.minCopyNumber() : -1,
                    cnData[CENTROMERE_CN],
                    arm == CHROMOSOME_ARM_P ? cnData[P_ARM_TELOMERE_CN] : cnData[Q_ARM_TELOMERE_CN]));

            writer.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing cluster-data to outputFile: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mFileWriter);
    }

    private void checkPseudoGeneAnnotations()
    {
        for(final SvCluster cluster : mClusters)
        {
            isSpecificCluster(cluster);

            for(final SvLinkedPair pair : cluster.getLinkedPairs())
            {
                if(pair.length() > SHORT_TI_LENGTH * 8)
                    continue;

                final SvBreakend lower = pair.getBreakend(true);
                final SvBreakend upper = pair.getBreakend(false);

                // for any TI falling within the same gene, check for an exon boundary match
                if(lower.getSV().getGenesList(lower.usesStart()).isEmpty() || upper.getSV().getGenesList(upper.usesStart()).isEmpty())
                    continue;

                for(final GeneAnnotation gene1 : lower.getSV().getGenesList(lower.usesStart()))
                {
                    for(final GeneAnnotation gene2 : upper.getSV().getGenesList(upper.usesStart()))
                    {
                        if(!gene1.GeneName.equals(gene2.GeneName))
                            continue;

                        final String exonData = mGeneTranscriptCollection.getExonDetailsForPosition(gene1, lower.position(), upper.position());

                        if(!exonData.isEmpty())
                        {
                            //LOGGER.info("sample({}) cluster({}) pair({}) matches gene({}) exon({})",
                            //        mSampleId, cluster.id(), pair.toString(), gene1.GeneName, exonData);

                            pair.setExonMatchData(exonData);
                        }
                    }
                }
            }
        }
    }
}
