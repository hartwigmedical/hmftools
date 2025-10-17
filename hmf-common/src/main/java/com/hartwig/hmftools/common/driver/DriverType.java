package com.hartwig.hmftools.common.driver;

import java.util.EnumSet;

public enum DriverType
{
    AMP,
    PARTIAL_AMP,
    DEL,
    HET_DEL,
    HOM_DUP_DISRUPTION,
    HOM_DEL_DISRUPTION,
    DISRUPTION,
    MUTATION,
    GERMLINE_MUTATION,
    GERMLINE_DELETION,
    GERMLINE_DISRUPTION,
    GERMLINE_HOM_DUP_DISRUPTION,
    UNKNOWN;

    public static final EnumSet<DriverType> DRIVERS_PURPLE_GERMLINE = EnumSet.of(GERMLINE_MUTATION, GERMLINE_DELETION);

    public static final EnumSet<DriverType> DRIVERS_PURPLE_SOMATIC = EnumSet.of(AMP, PARTIAL_AMP, DEL, HET_DEL, MUTATION);

    public static final EnumSet<DriverType> DRIVERS_LINX_GERMLINE = EnumSet.of(GERMLINE_DISRUPTION, GERMLINE_HOM_DUP_DISRUPTION);

    public static final EnumSet<DriverType> DRIVERS_LINX_SOMATIC = EnumSet.of(HOM_DUP_DISRUPTION, HOM_DEL_DISRUPTION, DISRUPTION);

    public static boolean isGermline(final DriverType type)
    {
        return type == GERMLINE_DELETION || type == GERMLINE_DISRUPTION || type == GERMLINE_MUTATION || type == GERMLINE_HOM_DUP_DISRUPTION;
    }

    public static DriverType checkConvertType(final String driverTypeStr)
    {
        // allow reading of old types for backwards compatibility
        try
        {
            return DriverType.valueOf(driverTypeStr);
        }
        catch(IllegalArgumentException e)
        {
            return UNKNOWN;
        }
    }

}
