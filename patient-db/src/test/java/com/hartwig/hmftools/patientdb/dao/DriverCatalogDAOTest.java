package com.hartwig.hmftools.patientdb.dao;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.drivercatalog.DriverType;

import org.junit.Test;

public class DriverCatalogDAOTest {

    @Test
    public void testDriverCatalogTypesAllSet() {

        for (DriverType value : DriverType.values()) {
            int count = 0;
            if (DriverCatalogDAO.DRIVER_CATALOG_GERMLINE.contains(value)) {
                count++;
            }
            if (DriverCatalogDAO.DRIVER_CATALOG_LINX.contains(value)) {
                count++;
            }
            if (DriverCatalogDAO.DRIVER_CATALOG_SOMATIC.contains(value)) {
                count++;
            }

            assertEquals(1, count);
        }
    }

}
