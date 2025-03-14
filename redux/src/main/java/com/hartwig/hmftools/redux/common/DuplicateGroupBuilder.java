package com.hartwig.hmftools.redux.common;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.umi.UmiConfig;
import com.hartwig.hmftools.redux.umi.UmiGroupBuilder;

import htsjdk.samtools.SAMRecord;

public class DuplicateGroupBuilder
{
    private final UmiConfig mUmiConfig;
    private final Statistics mStats;
    private final boolean mFormConsensus;
    private final UmiGroupBuilder mUmiGroupBuilder;

    public DuplicateGroupBuilder(final ReduxConfig config)
    {
        mFormConsensus = config.FormConsensus;
        mUmiConfig = config.UMIs;
        mStats = new Statistics();
        mUmiGroupBuilder = new UmiGroupBuilder(config.Sequencing, config.UMIs, mStats.UmiStats);
    }

    public Statistics statistics() { return mStats; }

    public static final Set<String> DBG_READ_NAMES = Sets.newHashSet(Lists.newArrayList(
            "A00260:691:H55JLDSXC:1:1243:10845:22686:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:1347:20672:36542:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:1360:2166:28855:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:2219:23185:8547:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:1309:6533:31594:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2152:24975:4820:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2464:11559:9236:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1152:3965:31501:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1318:23592:18145:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2163:26476:15499:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2216:7916:24846:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2375:25771:29434:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2505:26142:15515:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2531:11071:30780:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2627:22598:19398:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1520:10131:4946:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1461:4363:22138:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:2269:13675:16814:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:2338:17933:11256:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:2363:26784:5071:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:2603:30219:2284:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:2663:15637:14575:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:1619:11460:34084:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:1620:8377:3035:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2126:11803:6496:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2353:28691:2049:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2604:13783:15248:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1203:8323:25019:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1269:6451:6840:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1371:20320:34898:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2232:27995:36354:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2356:9471:23343:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2371:21730:27320:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1143:16721:10942:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1424:22119:5134:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1641:30725:15186:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:2306:21287:21292:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:2578:7419:28275:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:1130:14796:20603:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:1308:28293:10661:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:1407:11397:33160:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:2234:3070:25379:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:1424:13548:33285:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:1676:1542:33317:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2431:19922:22686:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2508:4752:25285:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2540:8196:25238:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2559:32904:24565:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1309:13919:17957:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1402:13259:23610:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1540:28917:19977:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1105:26096:14215:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1111:28402:28354:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1252:5909:21026:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1505:6153:18380:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1664:24858:21621:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1664:28782:23876:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:2556:31033:32471:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:1561:32533:7482:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1360:6696:3881:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:2359:8748:36620:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:2450:15528:20964:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:1243:10845:22686:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:1347:20672:36542:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:1360:2166:28855:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:1:2219:23185:8547:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:1309:6533:31594:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2152:24975:4820:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:2:2464:11559:9236:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1152:3965:31501:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:1318:23592:18145:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2163:26476:15499:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2216:7916:24846:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2375:25771:29434:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2505:26142:15515:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2531:11071:30780:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:3:2627:22598:19398:TAAGTGGC+ATCCGTTG",
    "A00260:691:H55JLDSXC:4:1520:10131:4946:TAAGTGGC+ATCCGTTG"
    ));

//    public static final Set<String> DBG_READ_NAMES = Set.of("A00260:691:H55JLDSXC:1:1243:10845:22686:TAAGTGGC+ATCCGTTG");

    public List<DuplicateGroup> processDuplicateGroups(
            final List<DuplicateGroup> rawDuplicateGroups, final List<ReadInfo> singleFragments, boolean captureStats)
    {
        if(mUmiConfig.Enabled)
        {
            List<DuplicateGroup> umiGroups = mUmiGroupBuilder.processUmiGroups(rawDuplicateGroups, singleFragments, captureStats);

            if(captureStats)
            {
                for(DuplicateGroup umiGroup: umiGroups)
                {
                    mStats.addUmiGroup(umiGroup.readCount(), umiGroup.hasDualStrand());
                }
            }

            List<DuplicateGroup> filteredUmiGroups = Lists.newArrayList();
            for(DuplicateGroup umiGroup : umiGroups)
            {
                boolean found = false;
                for(SAMRecord read : umiGroup.reads())
                {
                    if(DBG_READ_NAMES.contains(read.getReadName()))
                    {
                        found = true;
                        break;
                    }
                }

                if(found)
                {
                    filteredUmiGroups.add(umiGroup);
                }
            }

            List<ReadInfo> filteredSingleFragments = Lists.newArrayList();
            for(ReadInfo readInfo : singleFragments)
            {
                SAMRecord read = readInfo.read();
                if(DBG_READ_NAMES.contains(read.getReadName()))
                {
                    filteredSingleFragments.add(readInfo);
                }
            }

            singleFragments.clear();
            singleFragments.addAll(filteredSingleFragments);

            return filteredUmiGroups;
        }

        if(captureStats)
        {
            for(DuplicateGroup duplicateGroup : rawDuplicateGroups)
            {
                mStats.addDuplicateGroup(duplicateGroup.readCount());

                if(!mFormConsensus)
                    setPrimaryRead(duplicateGroup);
            }
        }

        return rawDuplicateGroups;
    }

    private static void setPrimaryRead(final DuplicateGroup duplicateGroup)
    {
        SAMRecord maxRead = null;
        double maxBaseQual = 0;

        for(SAMRecord read : duplicateGroup.reads())
        {
            double avgBaseQual = calcBaseQualAverage(read);

            if(avgBaseQual > maxBaseQual)
            {
                maxBaseQual = avgBaseQual;
                maxRead = read;
            }
        }

        duplicateGroup.setPrimaryRead(maxRead);
    }

    public static double calcBaseQualAverage(final SAMRecord read)
    {
        int readBaseCount = 0;
        int readBaseQualTotal = 0;

        for(int i = 0; i < read.getBaseQualities().length; ++i)
        {
            ++readBaseCount;
            readBaseQualTotal += read.getBaseQualities()[i];
        }

        return readBaseCount > 0 ? readBaseQualTotal / (double)readBaseCount : 0;
    }
}
