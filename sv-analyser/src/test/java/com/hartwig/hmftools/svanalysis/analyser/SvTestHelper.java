package com.hartwig.hmftools.svanalysis.analyser;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svanalysis.analysis.ClusterAnnotations.ALL_ANNOTATIONS;
import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.DEFAULT_PROXIMITY_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvSampleAnalyser.setSvCopyNumberData;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFactory;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.analysis.CNAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.ClusterAnalyser;
import com.hartwig.hmftools.svanalysis.analysis.LinkFinder;
import com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods;
import com.hartwig.hmftools.svanalysis.analysis.SvUtilities;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvaConfig;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;

public class SvTestHelper
{
    public String SampleId;
    public List<SvVarData> AllVariants;
    public SvaConfig Config;
    public SvClusteringMethods ClusteringMethods;
    public ClusterAnalyser Analyser;
    public CNAnalyser CopyNumberAnalyser;

    private int mNextVarId;

    public SvTestHelper()
    {
        Config = new SvaConfig(DEFAULT_PROXIMITY_DISTANCE);

        ClusteringMethods = new SvClusteringMethods(Config.ProximityDistance);
        Analyser = new ClusterAnalyser(Config, ClusteringMethods);
        CopyNumberAnalyser = new CNAnalyser("", null);
        Analyser.setCopyNumberAnalyser(CopyNumberAnalyser);

        Analyser.setRunValidationChecks(true);

        SampleId = "TEST";
        AllVariants = Lists.newArrayList();

        Analyser.setSampleData(SampleId, AllVariants);
        mNextVarId = 0;

        Configurator.setRootLevel(Level.DEBUG);
    }

    public final String nextVarId() { return String.format("%d", mNextVarId++); }
    public void logVerbose(boolean toggle)
    {
        Config.LogVerbose = toggle;
        Analyser.getChainFinder().setLogVerbose(toggle);
        Analyser.getLinkFinder().setLogVerbose(toggle);
    }

    public void preClusteringInit()
    {
        ClusteringMethods.populateChromosomeBreakendMap(AllVariants);

        addCopyNumberData();

        ClusteringMethods.annotateNearestSvData();
        LinkFinder.findDeletionBridges(ClusteringMethods.getChrBreakendMap());
        ClusteringMethods.setSimpleVariantLengths(SampleId);
    }

    public void addClusterAndSVs(final SvCluster cluster)
    {
        Analyser.getClusters().add(cluster);
        AllVariants.addAll(cluster.getSVs());
    }

    public void clearClustersAndSVs()
    {
        // in case SVs are to be used again and re-clustered
        for(SvVarData var : AllVariants)
        {
            var.setCluster(null);
        }

        AllVariants.clear();
        ClusteringMethods.clearLOHBreakendData(SampleId);
        Analyser.getClusters().clear();
    }

    public void mergeOnProximity()
    {
        ClusteringMethods.clusterByProximity(AllVariants, Analyser.getClusters());
    }


    public final List<SvCluster> getClusters() { return Analyser.getClusters(); }



    public static SvVarData createSv(final String varId, final String chrStart, final String chrEnd,
            long posStart, long posEnd, int orientStart, int orientEnd, StructuralVariantType type, final String insertSeq)
    {
        return createTestSv(varId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, type,
                2, 2, 1, 1, 1, insertSeq);
    }

    // for convenience
    public static SvVarData createDel(final String varId, final String chromosome, long posStart, long posEnd)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, 1, -1, DEL,
                2, 2, 1, 1, 1, "");
    }

    public static SvVarData createIns(final String varId, final String chromosome, long posStart, long posEnd)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, 1, -1, INS,
                2, 2, 1, 1, 1, "");
    }

    public static SvVarData createDup(final String varId, final String chromosome, long posStart, long posEnd)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, -1, 1, DUP,
                3, 3, 1, 1, 1, "");
    }

    public static SvVarData createInv(final String varId, final String chromosome, long posStart, long posEnd, int orientation)
    {
        return createTestSv(varId, chromosome, chromosome, posStart, posEnd, orientation, orientation, INV,
                orientation == 1 ? 4 : 3, orientation == 1 ? 3 : 4, 1, 1, 1, "");
    }

    public static SvVarData createSgl(final String varId, final String chromosome, long position, int orientation, boolean isNoneSegment)
    {
        SvVarData var = createTestSv(varId, chromosome, "0", position, -1, orientation, -1, SGL,
                3, 0, 1, 0, 1, "");

        return var;
    }

    public static SvVarData createBnd(final String varId, final String chrStart, long posStart, int orientStart, final String chrEnd, long posEnd, int orientEnd)
    {
        SvVarData var = createTestSv(varId, chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, BND,
                3, 3, 1, 1, 1, "");

        return var;
    }

    public static SvVarData createTestSv(final String varId, final String chrStart, final String chrEnd,
            long posStart, long posEnd, int orientStart, int orientEnd, StructuralVariantType type,
            double cnStart, double cnEnd, double cnChgStart, double cnChgEnd, double ploidy, final String insertSeq)
    {
        StructuralVariantData svData =
                ImmutableStructuralVariantData.builder()
                        .id(varId)
                        .startChromosome(chrStart)
                        .endChromosome(chrEnd)
                        .startPosition(posStart)
                        .endPosition(posEnd)
                        .startOrientation((byte)orientStart)
                        .endOrientation((byte)orientEnd)
                        .startAF(1.0)
                        .adjustedStartAF(1.0)
                        .adjustedStartCopyNumber(cnStart)
                        .adjustedStartCopyNumberChange(cnChgStart)
                        .endAF(1.0)
                        .adjustedEndAF(1.0)
                        .adjustedEndCopyNumber(cnEnd)
                        .adjustedEndCopyNumberChange(cnChgEnd)
                        .ploidy(ploidy)
                        .type(type)
                        .homology("")
                        .vcfId("")
                        .insertSequence(insertSeq)
                        .insertSequenceAlignments("")
                        .filter("PASS")
                        .imprecise(false)
                        .qualityScore(0.0)
                        .event("")
                        .startTumourVariantFragmentCount(0)
                        .startTumourReferenceFragmentCount(0)
                        .startNormalVariantFragmentCount(0)
                        .startNormalReferenceFragmentCount(0)
                        .endTumourVariantFragmentCount(0)
                        .endTumourReferenceFragmentCount(0)
                        .endNormalVariantFragmentCount(0)
                        .endNormalReferenceFragmentCount(0)
                        .startIntervalOffsetStart(0)
                        .startIntervalOffsetEnd(0)
                        .endIntervalOffsetStart(0)
                        .endIntervalOffsetEnd(0)
                        .inexactHomologyOffsetStart(0)
                        .inexactHomologyOffsetEnd(0)
                        .startLinkedBy("")
                        .endLinkedBy("")
                        .startRefContext("")
                        .endRefContext("")
                        .recovered(false)
                        .insertSequenceRepeatClass("")
                        .insertSequenceRepeatType("")
                        .insertSequenceRepeatOrientation((byte)0)
                        .insertSequenceRepeatCoverage(0.0)
                        .build();

        SvVarData var = new SvVarData(svData);

        String startArm = getChromosomalArm(var.chromosome(true), var.position(true));

        String endArm;
        if(!var.isNullBreakend())
            endArm = getChromosomalArm(var.chromosome(false), var.position(false));
        else
            endArm = CHROMOSOME_ARM_P;

        var.setChromosomalArms(startArm, endArm);

        // by default
        var.setPloidyRecalcData(var.getSvData().ploidy(), var.getSvData().ploidy());

        return var;
    }

    public void setDefaultPloidyCalcData()
    {
        AllVariants.stream().forEach(x -> x.setPloidyRecalcData(x.getSvData().ploidy(), x.getSvData().ploidy()));
    }

    public void addCopyNumberData()
    {
        // use SV breakend data to re-create the copy number segments
        // NOTE: positions adjusted for orientation are not done correctly
        Map<String, List<SvCNData>> chrCnDataMap = CopyNumberAnalyser.getChrCnDataMap();
        final Map<String, List<SvBreakend>> chrBreakendMap = ClusteringMethods.getChrBreakendMap();
        Map<String,SvCNData[]> svIdCnDataMap = CopyNumberAnalyser.getSvIdCnDataMap();

        chrCnDataMap.clear();
        svIdCnDataMap.clear();

        double nonDisruptedAP = 1;

        int cnId = 0;
        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final String chromosome = entry.getKey();
            List<SvBreakend> breakendList = entry.getValue();
            List<SvCNData> cnDataList = Lists.newArrayList();
            chrCnDataMap.put(chromosome, cnDataList);

            long centromerePosition = SvUtilities.getChromosomalArmLength(chromosome, CHROMOSOME_ARM_P);
            long chromosomeLength = SvUtilities.CHROMOSOME_LENGTHS.get(chromosome);

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvVarData var = breakend.getSV();

                SvCNData cnData = null;

                if (i == 0)
                {
                    // add telomere segment at start, and centromere as soon as the breakend crosses the centromere
                    if(breakend.arm() == CHROMOSOME_ARM_Q)
                    {
                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, 0, centromerePosition,
                                breakend.getCopyNumber(true),
                                TELOMERE.toString(), CENTROMERE.toString(),
                                1, 0.5, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);

                        extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, breakend.position() - 1,
                                breakend.getCopyNumber(true),
                                CENTROMERE.toString(), var.type().toString(),
                                1, 0.5, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, 0, breakend.position() - 1,
                                breakend.getCopyNumber(true),
                                TELOMERE.toString(), var.type().toString(),
                                1, 0.5, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                }

                if (i < breakendList.size() - 1)
                {
                    final SvBreakend nextBreakend = breakendList.get(i + 1);

                    if(breakend.arm() == CHROMOSOME_ARM_P && nextBreakend.arm() == CHROMOSOME_ARM_Q)
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), centromerePosition-1,
                                breakend.getCopyNumber(false),
                                var.type().toString(), CENTROMERE.toString(),
                                1, 0.5, 100);

                        cnData.setIndex(cnDataList.size());
                        cnDataList.add(cnData);

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, nextBreakend.position() - 1,
                                breakend.getCopyNumber(false),
                                CENTROMERE.toString(), nextBreakend.getSV().type().toString(),
                                1, 0.5, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), nextBreakend.position() - 1,
                                breakend.getCopyNumber(false),
                                var.type().toString(), nextBreakend.getSV().type().toString(),
                                1, 0.5, 100);

                        cnData.setIndex(cnDataList.size());
                        cnDataList.add(cnData);
                    }
                }
                else
                {
                    // last breakend runs out to the telomere
                    if(breakend.arm() == CHROMOSOME_ARM_P)
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), centromerePosition - 1,
                                breakend.getCopyNumber(false),
                                var.type().toString(), CENTROMERE.toString(),
                                1, 0.5, 100);

                        cnData.setIndex(cnDataList.size());
                        cnDataList.add(cnData);

                        SvCNData extraCnData = new SvCNData(cnId++, chromosome, centromerePosition, chromosomeLength,
                                breakend.getCopyNumber(false),
                                CENTROMERE.toString(), TELOMERE.toString(),
                                1, 0.5, 100);

                        extraCnData.setIndex(cnDataList.size());
                        cnDataList.add(extraCnData);
                    }
                    else
                    {
                        cnData = new SvCNData(cnId++, chromosome, breakend.position(), chromosomeLength,
                                breakend.getCopyNumber(false),
                                var.type().toString(), TELOMERE.toString(),
                                1, 0.5, 100);

                        cnData.setIndex(cnDataList.size());
                        cnDataList.add(cnData);
                    }
                }

                SvCNData[] cnDataPair = svIdCnDataMap.get(var.id());

                if(cnDataPair == null)
                {
                    cnDataPair = new SvCNData[2];
                    svIdCnDataMap.put(var.id(), cnDataPair);
                }

                cnDataPair[breakend.usesStart() ? SVI_START : SVI_END] = cnData;
            }
        }

        setSvCopyNumberData(
                AllVariants,
                CopyNumberAnalyser.getSampleSvPloidyCalcMap().get(SampleId),
                CopyNumberAnalyser.getSvIdCnDataMap(),
                CopyNumberAnalyser.getChrCnDataMap());

    }

}
