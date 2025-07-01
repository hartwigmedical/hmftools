package com.hartwig.hmftools.chord.prep;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.indel.IndelPrep;
import com.hartwig.hmftools.chord.snv.SnvPrep;
import com.hartwig.hmftools.chord.sv.SvPrep;

import org.apache.logging.log4j.Level;
import org.jetbrains.annotations.Nullable;

public class SamplePrepTask implements Callable<Void>
{
    private final ChordConfig mConfig;

    private final int mSampleIndex;
    private final String mSampleId;

    public List<MutContextCount> mContextCounts = new ArrayList<>();

    @Nullable
    private final ConcurrentHashMap<String, List<MutContextCount>> mContextCountsMatrix;

    private static final int FEW_SAMPLES_THRESHOLD = 10;
    private static final int PROGRESS_INTERVAL = 1000;

    public SamplePrepTask(
            final ChordConfig chordConfig,
            final int sampleIndex,
            @Nullable final ConcurrentHashMap<String, List<MutContextCount>> contextCountsMatrix
    )
    {
        mConfig = chordConfig;

        mSampleIndex = sampleIndex;
        mSampleId = mConfig.SampleIds.get(sampleIndex);

        mContextCountsMatrix = contextCountsMatrix;
    }

    public void processSample()
    {
        String logPrefix = mConfig.isMultiSample() ?
                String.format("%s/%s: sample(%s): ", mSampleIndex + 1, mConfig.SampleIds.size(), mSampleId) :
                "";

        boolean showCoarseProgress = mConfig.isMultiSample() &&
                CHORD_LOGGER.getLevel().isMoreSpecificThan(Level.INFO) &&
                (mSampleIndex < FEW_SAMPLES_THRESHOLD || mSampleIndex % PROGRESS_INTERVAL == 0);

        if(showCoarseProgress)
        {
            CHORD_LOGGER.info(logPrefix + "Running " + ChordDataPrep.class.getSimpleName());
        }

        List<VariantTypePrep<?>> variantTypePreps = List.of(
                new SnvPrep(mConfig).logPrefix(logPrefix),
                new IndelPrep(mConfig).logPrefix(logPrefix),
                new SvPrep(mConfig).logPrefix(logPrefix)
        );

        for(VariantTypePrep variantTypePrep : variantTypePreps)
        {
            List<MutContextCount> variantTypeContextCounts = variantTypePrep.countMutationContexts(mSampleId);

            mContextCounts.addAll(variantTypeContextCounts);
        }

        if(mConfig.isMultiSample())
        {
            synchronized (mContextCountsMatrix)
            {
                mContextCountsMatrix.put(mSampleId, mContextCounts);
            }

            mContextCounts = null;
        }
    }

    @Override
    public Void call()
    {
        processSample();
        return null;
    }
}
