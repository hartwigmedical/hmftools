package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicUser {

    @NotNull
    public abstract String username();

    @NotNull
    public abstract String areaOfExpertise();

    @NotNull
    public abstract CivicOrganization organization();

    @NotNull
    public abstract String twitterHandle();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String bio();

    @NotNull
    public abstract String url();

    @NotNull
    public abstract String createdAt();

    @NotNull
    public abstract CivicAvatars avatars();

    @NotNull
    public abstract String acceptedLicense();

    @NotNull
    public abstract String affiliation();

    @NotNull
    public abstract String avatarUrl();

    @NotNull
    public abstract String role();

    @NotNull
    public abstract String facebookProfile();

    @NotNull
    public abstract String linkedinProfile();

    @NotNull
    public abstract String orcid();

    @NotNull
    public abstract String displayName();

    @NotNull
    public abstract String lastSeenAt();

    @NotNull
    public abstract String featuredExpert();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String signupComplete();
}
