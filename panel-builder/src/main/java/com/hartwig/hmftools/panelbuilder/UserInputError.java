package com.hartwig.hmftools.panelbuilder;

// Fatal error which is directly caused by user input and must be rectified by the user.
public class UserInputError extends RuntimeException
{
    public UserInputError(String message)
    {
        super(message);
    }
}
