package com.hartwig.hmftools.peach.data_loader;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.peach.effect.DrugInfo;
import com.hartwig.hmftools.peach.effect.ImmutableDrugInfo;

import org.junit.Test;

public class DrugInfoLoaderTest
{
    @Test
    public void testLoad()
    {
        String filePath = getTestResourcePath("drugs.tsv");
        List<DrugInfo> drugInfos = DrugInfoLoader.loadDrugInfos(filePath);

        assertEquals(4, drugInfos.size());

        DrugInfo expectedDrugInfo0 = ImmutableDrugInfo.builder()
            .drugName("5-Fluorouracil")
            .geneName("DPYD")
            .prescriptionInfoUrl("https://www.pharmgkb.org/guidelineAnnotation/PA166104939")
            .build();
        assertEquals(expectedDrugInfo0, drugInfos.get(0));

        DrugInfo expectedDrugInfo1 = ImmutableDrugInfo.builder()
                .drugName("Capecitabine")
                .geneName("DPYD")
                .prescriptionInfoUrl("https://www.pharmgkb.org/guidelineAnnotation/PA166104963")
                .build();
        assertEquals(expectedDrugInfo1, drugInfos.get(1));

        DrugInfo expectedDrugInfo2 = ImmutableDrugInfo.builder()
                .drugName("Tegafur")
                .geneName("DPYD")
                .prescriptionInfoUrl("https://www.pharmgkb.org/guidelineAnnotation/PA166104944")
                .build();
        assertEquals(expectedDrugInfo2, drugInfos.get(2));

        DrugInfo expectedDrugInfo3 = ImmutableDrugInfo.builder()
                .drugName("Irinotecan")
                .geneName("UGT1A1")
                .prescriptionInfoUrl("https://www.pharmgkb.org/guidelineAnnotation/PA166104951")
                .build();
        assertEquals(expectedDrugInfo3, drugInfos.get(3));
    }
}
