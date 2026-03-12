package com.hartwig.hmftools.qsee.prep.category.msindel;

public enum RepeatType
{
    REPEAT_A_T("A/T repeat"),
    REPEAT_C_G("C/G repeat"),
    REPEAT_2_BP("2bp repeat"),
    REPEAT_3_BP_OR_LONGER(">=3bp repeat");

    private final String mDisplayName;

    RepeatType(String displayName)
    {
        mDisplayName = displayName;
    }

    public String displayName() { return mDisplayName; }

    public static RepeatType fromDisplayName(String displayName)
    {
        if(displayName.equals(REPEAT_A_T.displayName()))
            return REPEAT_A_T;

        if(displayName.equals(REPEAT_C_G.displayName()))
            return REPEAT_C_G;

        if(displayName.equals(REPEAT_2_BP.displayName()))
            return REPEAT_2_BP;

        if(displayName.equals(REPEAT_3_BP_OR_LONGER.displayName()))
            return REPEAT_3_BP_OR_LONGER;

        throw new IllegalArgumentException("Unknown repeat type: " + displayName);
    }

    public static RepeatType fromRepeatUnit(String repeatUnit)
    {
        if(repeatUnit.equals("A/T"))
            return REPEAT_A_T;

        if(repeatUnit.equals("C/G"))
            return REPEAT_C_G;

        else if(repeatUnit.matches("^\\w{2}/.*"))
            return RepeatType.REPEAT_2_BP;

        else if(repeatUnit.matches("^\\d+bp repeat"))
            return RepeatType.REPEAT_3_BP_OR_LONGER;

        else
            throw new IllegalArgumentException("Unexpected repeat unit: " + repeatUnit);
    }
}
