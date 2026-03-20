package com.hartwig.hmftools.finding;

import java.io.IOException;

import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.datamodel.TestOrangeJsonWriter;

import org.junit.Ignore;
import org.junit.Test;

public class FindingFactoryTest
{
    @Ignore
    @Test
    public void fromOrangeRecord() throws IOException
    {
        FindingRecordFactory.fromOrangeRecord(TestOrangeJsonWriter.createOrangeRecord(), null, null, Gender.FEMALE);
    }
}
