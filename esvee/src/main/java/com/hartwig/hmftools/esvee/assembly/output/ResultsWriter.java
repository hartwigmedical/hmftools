package com.hartwig.hmftools.esvee.assembly.output;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;

import java.io.BufferedWriter;

import com.hartwig.hmftools.esvee.AssemblyConfig;
import com.hartwig.hmftools.esvee.alignment.DecoyChecker;
import com.hartwig.hmftools.esvee.utils.TruthsetAnnotation;

public class ResultsWriter
{
    private final BufferedWriter mDecoyMatchWriter;
    private final AssemblyWriter mAssemblyWriter;
    private final BreakendWriter mBreakendWriter;
    private final AssemblyReadWriter mReadWriter;
    private final PhaseGroupBuildWriter mPhaseGroupBuildWriter;
    private final BamWriter mBamWriter;
    private final TruthsetAnnotation mTruthsetAnnotation;

    public ResultsWriter(final AssemblyConfig config)
    {
        mTruthsetAnnotation = new TruthsetAnnotation(config.TruthsetFile);
        mAssemblyWriter = new AssemblyWriter(config, mTruthsetAnnotation);
        mBreakendWriter = new BreakendWriter(config, mTruthsetAnnotation);
        mReadWriter = new AssemblyReadWriter(config);
        mPhaseGroupBuildWriter = new PhaseGroupBuildWriter(config);
        mBamWriter = new BamWriter(config);
        mDecoyMatchWriter = DecoyChecker.initialiseWriter(config);
    }

    public BufferedWriter decoyMatchWriter() { return mDecoyMatchWriter; }
    public AssemblyWriter assemblyWriter() { return mAssemblyWriter; }
    public BreakendWriter breakendWriter() { return mBreakendWriter; }
    public AssemblyReadWriter readWriter() { return mReadWriter; }
    public PhaseGroupBuildWriter phaseGroupBuildWriter() { return mPhaseGroupBuildWriter; }
    public BamWriter bamWriter() { return mBamWriter; }

    public void close()
    {
        mAssemblyWriter.close();
        mBreakendWriter.close();
        mReadWriter.close();
        mPhaseGroupBuildWriter.close();
        mBamWriter.close();
        closeBufferedWriter(mDecoyMatchWriter);
    }
}
