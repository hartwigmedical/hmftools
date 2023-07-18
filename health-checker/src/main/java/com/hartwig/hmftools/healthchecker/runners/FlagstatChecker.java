package com.hartwig.hmftools.healthchecker.runners;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.flagstat.FlagstatFile;
import com.hartwig.hmftools.healthchecker.result.QCValue;
import com.hartwig.hmftools.healthchecker.result.QCValueType;

public class FlagstatChecker extends FileBasedHealthChecker<Flagstat>
{
    public FlagstatChecker(
            final String flagstatFile, final QCValueType proportionMapped, final QCValueType proportionDuplicate) throws IOException
    {
        super(FlagstatFile.read(flagstatFile),
                flagstat -> List.of(
                        new QCValue(proportionMapped, String.valueOf(flagstat.mappedProportion())),
                        new QCValue(proportionDuplicate, String.valueOf(flagstat.duplicateProportion()))));

    }
}
