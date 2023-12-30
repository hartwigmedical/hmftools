package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.SvConfig.SV_LOGGER;

import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.SvConfig;
import com.hartwig.hmftools.esvee.common.JunctionGroup;
import com.hartwig.hmftools.esvee.processor.PrimaryAssemblyResult;
import com.hartwig.hmftools.esvee.read.BamReader;

public class AssemblyThread extends Thread
{
    private final SvConfig mConfig;
    private final Queue<JunctionGroup> mJunctionGroups;
    private final int mJunctionCount;

    private final BamReader mBamReader;

    private final List<PrimaryAssemblyResult> mPrimaryAssemblyResults;

    public AssemblyThread(final SvConfig config, final Queue<JunctionGroup> junctionGroups)
    {
        mConfig = config;
        mJunctionGroups = junctionGroups;
        mJunctionCount = junctionGroups.size();

        mBamReader = new BamReader(config);

        mPrimaryAssemblyResults = Lists.newArrayList();

        start();
    }

    public List<PrimaryAssemblyResult> primaryAssemblyResults() { return mPrimaryAssemblyResults; }

    @Override
    public void run()
    {
        while(true)
        {
            try
            {
                int remainingCount = mJunctionGroups.size();
                int processedCount = mJunctionCount - remainingCount;

                JunctionGroup junctionGroup = mJunctionGroups.remove();

                if(processedCount > 0 && (processedCount % 100) == 0)
                {
                    SV_LOGGER.info("processed {} junction groups, remaining({})", processedCount, remainingCount);
                }

                JunctionGroupAssembler primaryAssembler = new JunctionGroupAssembler(mConfig, mBamReader, junctionGroup);
                primaryAssembler.run();

                mPrimaryAssemblyResults.addAll(primaryAssembler.primaryAssemblyResults());
            }
            catch(NoSuchElementException e)
            {
                SV_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }
}
