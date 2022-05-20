package com.hartwig.hmftools.patientreporter.remarks;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class SpecialRemarkFileTest {

    private static final String SAMPLE_SPECIAL_REMARK_TSV = Resources.getResource("special_remark/sample_special_remark.tsv").getPath();

    @Test
    public void specialRemarkFromTsvWithNewLines() throws IOException {
        SpecialRemarkModel specialRemarkModel = SpecialRemarkFile.buildFromTsv(SAMPLE_SPECIAL_REMARK_TSV);
        assertEquals(1, specialRemarkModel.specialRemarkCount());

        String specialRemark = specialRemarkModel.findSpecialRemarkForSample("sample");

        assertEquals(3, specialRemark.split("\n").length);
    }
}