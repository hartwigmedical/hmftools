package com.hartwig.hmftools.patientdb.dao;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class StructuralVariantDaoTest
{

    @Test
    public void testLimitSizeOfCSV()
    {
        String one = "asm12/23";
        String two = one + ",asm13/24";
        String three = two + ",asm12233/1";

        assertEquals(three, StructuralVariantDAO.limitSizeOfCSV(three.length(), three));

        for(int i = three.length() - 1; i >= two.length(); i--)
        {
            assertEquals(two, StructuralVariantDAO.limitSizeOfCSV(i, three));
        }

        for(int i = two.length() - 1; i >= one.length(); i--)
        {
            assertEquals(one, StructuralVariantDAO.limitSizeOfCSV(i, three));
        }
    }
}
