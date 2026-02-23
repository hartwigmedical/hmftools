package com.hartwig.hmftools.qsee.status;

public enum ComparisonOperator
{
    LESS_THAN("<"),
    LESS_THAN_OR_EQUAL("≤"),
    GREATER_THAN_OR_EQUAL("≥"),
    GREATER_THAN(">");

    public final String mOperatorString;

    ComparisonOperator(String operatorString)
    {
        mOperatorString = operatorString;
    }

    public String operatorString() { return mOperatorString; }
}
