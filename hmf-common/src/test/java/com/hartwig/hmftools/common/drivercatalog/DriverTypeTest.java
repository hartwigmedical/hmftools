package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class DriverTypeTest {

    @Test
    public void testDriverCatalogTypesAllSet() {
        for (DriverType value : DriverType.values()) {
            int count = 0;
            if (DriverType.DRIVERS_COPY_NUMBER.contains(value)) {
                count++;
            }
            if (DriverType.DRIVERS_GERMLINE.contains(value)) {
                count++;
            }
            if (DriverType.DRIVERS_LINX.contains(value)) {
                count++;
            }
            if (DriverType.DRIVERS_MUTATION.contains(value)) {
                count++;
            }

            assertEquals(1, count);
        }
    }
}
