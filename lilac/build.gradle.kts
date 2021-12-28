
plugins {
    id("com.hartwig.java-conventions")
    id("com.hartwig.build-version")
}

description = "HMF Tools - Lilac"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    testImplementation(libs.junit)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.lilac.LilacApplication")
}
