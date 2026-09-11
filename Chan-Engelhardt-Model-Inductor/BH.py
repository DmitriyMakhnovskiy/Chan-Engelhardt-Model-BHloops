
from Parameters import parameters as par

pi = 3.1415926535897932384626433832795
mu0 = 4 * pi * 1.0e-7  # permeability of vacuum (H/m)

Bs = par["Bs"]
Br = par["Br"]
Hc = par["Hc"]

# ============================================================
# MAJOR LOOP AND AUXILIARY FUNCTIONS
# ============================================================

def Upper(H):
    return Bs * (H + Hc) / (
        abs(H + Hc) + Hc * (Bs / Br - 1.0)
    ) + mu0 * H


def Lower(H):
    return Bs * (H - Hc) / (
        abs(H - Hc) + Hc * (Bs / Br - 1.0)
    ) + mu0 * H


def Binit(H):
    return 0.5 * (Upper(H) + Lower(H))


def UpperBmag(H):
    return Lower(H) + Br if H < 0.0 else Binit(H)


def LowerBmag(H):
    return Upper(H) - Br if H > 0.0 else Binit(H)


def BlatchUpper(H, scale, offset):
    return scale * UpperBmag(H) + offset


def BlatchLower(H, scale, offset):
    return scale * LowerBmag(H) + offset


def BlatchUpperLimit(H, Bd_turn):
    return Upper(H) - Bd_turn


def BlatchLowerLimit(H, Bd_turn):
    return Lower(H) + Bd_turn


# ============================================================
# FIG. 3A
# STEP 312.2
# SOLVE FOR H0 FOR INCREASING H
# ============================================================

def solve_H0_upper(Hlast, Bd_turn, tol=1e-10, max_iter=100):
    """
    Step 312.2 of Fig. 3A.

    Find H0 from:

        Binit(H0) = BlatchUpperLimit(H0)

    where:

        BlatchUpperLimit(H) = Upper(H) - Bd_turn

    The root is searched to the right of Hlast because
    this branch corresponds to increasing H.
    """

    def f(H):
        return Binit(H) - BlatchUpperLimit(H, Bd_turn)

    H_left = Hlast
    f_left = f(H_left)

    if abs(f_left) < tol:
        return H_left

    # Find a right boundary that brackets the root.
    step = max(1.0, 0.1 * abs(Hlast))

    H_right = H_left + step
    f_right = f(H_right)

    for _ in range(max_iter):

        if f_left * f_right <= 0.0:
            break

        step *= 2.0

        H_right = H_left + step
        f_right = f(H_right)

    else:
        raise RuntimeError(
            "Unable to bracket the root H0"
        )

    # Bisection method.
    for _ in range(max_iter):

        H_mid = 0.5 * (H_left + H_right)
        f_mid = f(H_mid)

        if abs(f_mid) < tol:
            return H_mid

        if f_left * f_mid <= 0.0:

            H_right = H_mid
            f_right = f_mid

        else:

            H_left = H_mid
            f_left = f_mid

    return 0.5 * (H_left + H_right)


# ============================================================
# FIG. 3A
# STEP 312.3
# SOLVE AFFINE MAP CONSTANTS FOR INCREASING H
# ============================================================

def solve_affine_upper(Hlast, Blast, H0):
    """
    Step 312.3 of Fig. 3A.

    The affine map is:

        BlatchUpper(H) =
            scale * UpperBmag(H) + offset

    The transformed curve passes through:

        (Hlast, Blast)

    and:

        (H0, Binit(H0))
    """

    x1 = UpperBmag(Hlast)
    y1 = Blast

    x2 = UpperBmag(H0)
    y2 = Binit(H0)

    denominator = x2 - x1

    if abs(denominator) < 1e-15:
        raise RuntimeError(
            "Unable to determine affine-map constants: "
            "UpperBmag(H0) equals UpperBmag(Hlast)"
        )

    scale = (y2 - y1) / denominator
    offset = y1 - scale * x1

    return scale, offset


# ============================================================
# FIG. 3B
# STEP 352.2
# SOLVE FOR H0 FOR DECREASING H
# ============================================================

def solve_H0_lower(Hlast, Bd_turn, tol=1e-10, max_iter=100):
    """
    Step 352.2 of Fig. 3B.

    Find H0 from:

        Binit(H0) = BlatchLowerLimit(H0)

    where:

        BlatchLowerLimit(H) = Lower(H) + Bd_turn

    The root is searched to the left of Hlast because
    this branch corresponds to decreasing H.
    """

    def f(H):
        return Binit(H) - BlatchLowerLimit(H, Bd_turn)

    H_right = Hlast
    f_right = f(H_right)

    if abs(f_right) < tol:
        return H_right

    # Find a left boundary that brackets the root.
    step = max(1.0, 0.1 * abs(Hlast))

    H_left = H_right - step
    f_left = f(H_left)

    for _ in range(max_iter):

        if f_left * f_right <= 0.0:
            break

        step *= 2.0

        H_left = H_right - step
        f_left = f(H_left)

    else:
        raise RuntimeError(
            "Unable to bracket the root H0"
        )

    # Bisection method.
    for _ in range(max_iter):

        H_mid = 0.5 * (H_left + H_right)
        f_mid = f(H_mid)

        if abs(f_mid) < tol:
            return H_mid

        if f_left * f_mid <= 0.0:

            H_right = H_mid
            f_right = f_mid

        else:

            H_left = H_mid
            f_left = f_mid

    return 0.5 * (H_left + H_right)


# ============================================================
# FIG. 3B
# STEP 352.3
# SOLVE AFFINE MAP CONSTANTS FOR DECREASING H
# ============================================================

def solve_affine_lower(Hlast, Blast, H0):
    """
    Step 352.3 of Fig. 3B.

    The affine map is:

        BlatchLower(H) =
            scale * LowerBmag(H) + offset

    The transformed curve passes through:

        (Hlast, Blast)

    and:

        (H0, Binit(H0))
    """

    x1 = LowerBmag(Hlast)
    y1 = Blast

    x2 = LowerBmag(H0)
    y2 = Binit(H0)

    denominator = x2 - x1

    if abs(denominator) < 1e-15:
        raise RuntimeError(
            "Unable to determine affine-map constants: "
            "LowerBmag(H0) equals LowerBmag(Hlast)"
        )

    scale = (y2 - y1) / denominator
    offset = y1 - scale * x1

    return scale, offset


# ============================================================
# BH_increasing
# CURRENT H IS INCREASING
# FIG. 3A
# ============================================================

def BH_increasing(
    H,
    was_inc,
    Hlast,
    Blast,
    OnInitMag,
    IsLatched,
    Bd,
    Bd_turn,
    H0,
    scale,
    offset
):

    # --------------------------------------------------------
    # Preserve all state variables.
    # --------------------------------------------------------

    OnInitMag_new = OnInitMag
    IsLatched_new = IsLatched

    Bd_new = Bd
    Bd_turn_new = Bd_turn

    H0_new = H0
    scale_new = scale
    offset_new = offset


    # ========================================================
    # STEP 308
    # WAS H INCREASING?
    # ========================================================

    if not was_inc:

        # H was decreasing and is now increasing.
        # A reversal occurred at (Hlast, Blast).

        # ====================================================
        # STEP 310
        #
        # IS Hlast <= 0
        # AND
        # Blast - Lower(Hlast) <= Br ?
        # ====================================================

        if (
            Hlast <= 0.0
            and Blast - Lower(Hlast) <= Br
        ):

            # ------------------------------------------------
            # STEP 314
            # CLEAR ISLATCHED
            # ------------------------------------------------

            IsLatched_new = False

        else:

            # ------------------------------------------------
            # STEP 312
            # CONSTRUCT THE AFFINE-MAP SOLUTION
            # ------------------------------------------------

            # STEP 312.1
            # SET Bd_turn = Bd

            if Bd is None:
                raise RuntimeError(
                    "Step 312 requires the previous Bd, "
                    "but Bd is undefined"
                )

            Bd_turn_new = Bd


            # STEP 312.2
            # SOLVE FOR H0

            H0_new = solve_H0_upper(
                Hlast,
                Bd_turn_new
            )


            # STEP 312.3
            # SOLVE AFFINE MAP CONSTANTS

            scale_new, offset_new = solve_affine_upper(
                Hlast,
                Blast,
                H0_new
            )


            # STEP 312.4
            # SET ISLATCHED

            IsLatched_new = True


        # ====================================================
        # STEP 316
        #
        # SET Bd = Blast - Lower(Hlast)
        # CLEAR ONINITMAG
        # ====================================================

        Bd_new = Blast - Lower(Hlast)

        OnInitMag_new = False


    # ========================================================
    # STEP 318
    # IS ONINITMAG SET?
    # ========================================================

    if OnInitMag_new:

        # STEP 332
        # B = Binit(H)

        B = Binit(H)

    else:

        # ====================================================
        # STEP 320
        #
        # IS ISLATCHED SET
        # AND
        # IS H < H0 ?
        # ====================================================

        use_latched_branch = False

        if IsLatched_new:

            if H0_new is None:
                raise RuntimeError(
                    "IsLatched is True, but H0 is undefined"
                )

            if H < H0_new:
                use_latched_branch = True


        if use_latched_branch:

            # ------------------------------------------------
            # STEP 324
            #
            # B = MIN(
            #     BlatchUpper(H),
            #     BlatchUpperLimit(H)
            # )
            # ------------------------------------------------

            if scale_new is None or offset_new is None:
                raise RuntimeError(
                    "IsLatched is True, "
                    "but scale or offset is undefined"
                )

            if Bd_turn_new is None:
                raise RuntimeError(
                    "IsLatched is True, "
                    "but Bd_turn is undefined"
                )

            B = min(
                BlatchUpper(
                    H,
                    scale_new,
                    offset_new
                ),
                BlatchUpperLimit(
                    H,
                    Bd_turn_new
                )
            )

        else:

            # ------------------------------------------------
            # STEP 322
            # CLEAR ISLATCHED
            # ------------------------------------------------

            IsLatched_new = False


            if Bd_new is None:
                raise RuntimeError(
                    "Step 326 requires Bd, "
                    "but Bd is undefined"
                )


            # =================================================
            # STEP 326
            #
            # IS Lower(H) + Bd > Binit(H)
            # AND
            # H >= 0 ?
            # =================================================

            if (
                Lower(H) + Bd_new > Binit(H)
                and H >= 0.0
            ):

                # STEP 328
                # SET ONINITMAG

                OnInitMag_new = True

                # STEP 332
                # B = Binit(H)

                B = Binit(H)

            else:

                # STEP 330
                # B = Lower(H) + Bd

                B = Lower(H) + Bd_new


    # ========================================================
    # RETURN UPDATED STATE
    # ========================================================

    return (
        B,
        OnInitMag_new,
        IsLatched_new,
        Bd_new,
        Bd_turn_new,
        H0_new,
        scale_new,
        offset_new
    )


# ============================================================
# BH_decreasing
# CURRENT H IS DECREASING
# FIG. 3B
# ============================================================

def BH_decreasing(
    H,
    was_dec,
    Hlast,
    Blast,
    OnInitMag,
    IsLatched,
    Bd,
    Bd_turn,
    H0,
    scale,
    offset
):

    # --------------------------------------------------------
    # Preserve all state variables.
    # --------------------------------------------------------

    OnInitMag_new = OnInitMag
    IsLatched_new = IsLatched

    Bd_new = Bd
    Bd_turn_new = Bd_turn

    H0_new = H0
    scale_new = scale
    offset_new = offset


    # ========================================================
    # STEP 348
    # WAS H DECREASING?
    # ========================================================

    if not was_dec:

        # H was increasing and is now decreasing.
        # A reversal occurred at (Hlast, Blast).

        # ====================================================
        # STEP 350
        #
        # IS Hlast >= 0
        # AND
        # Upper(Hlast) - Blast <= Br ?
        # ====================================================

        if (
            Hlast >= 0.0
            and Upper(Hlast) - Blast <= Br
        ):

            # ------------------------------------------------
            # STEP 354
            # CLEAR ISLATCHED
            # ------------------------------------------------

            IsLatched_new = False

        else:

            # ------------------------------------------------
            # STEP 352
            # CONSTRUCT THE AFFINE-MAP SOLUTION
            # ------------------------------------------------

            # STEP 352.1
            # SET Bd_turn = Bd

            if Bd is None:
                raise RuntimeError(
                    "Step 352 requires the previous Bd, "
                    "but Bd is undefined"
                )

            Bd_turn_new = Bd


            # STEP 352.2
            #
            # Binit(H0) = BlatchLowerLimit(H0)
            #

            H0_new = solve_H0_lower(
                Hlast,
                Bd_turn_new
            )


            # STEP 352.3
            # SOLVE AFFINE MAP CONSTANTS

            scale_new, offset_new = solve_affine_lower(
                Hlast,
                Blast,
                H0_new
            )


            # STEP 352.4
            # SET ISLATCHED

            IsLatched_new = True


        # ====================================================
        # STEP 356
        #
        # SET Bd = Upper(Hlast) - Blast
        # CLEAR ONINITMAG
        # ====================================================

        Bd_new = Upper(Hlast) - Blast

        OnInitMag_new = False


    # ========================================================
    # STEP 358
    # IS ONINITMAG SET?
    # ========================================================

    if OnInitMag_new:

        # STEP 372
        # B = Binit(H)

        B = Binit(H)

    else:

        # ====================================================
        # STEP 360
        #
        # IS ISLATCHED SET
        # AND
        # IS H > H0 ?
        # ====================================================

        use_latched_branch = False

        if IsLatched_new:

            if H0_new is None:
                raise RuntimeError(
                    "IsLatched is True, but H0 is undefined"
                )

            if H > H0_new:
                use_latched_branch = True


        if use_latched_branch:

            # ------------------------------------------------
            # STEP 364
            #
            # B = MAX(
            #     BlatchLower(H),
            #     BlatchLowerLimit(H)
            # )
            # ------------------------------------------------

            if scale_new is None or offset_new is None:
                raise RuntimeError(
                    "IsLatched is True, "
                    "but scale or offset is undefined"
                )

            if Bd_turn_new is None:
                raise RuntimeError(
                    "IsLatched is True, "
                    "but Bd_turn is undefined"
                )

            B = max(
                BlatchLower(
                    H,
                    scale_new,
                    offset_new
                ),
                BlatchLowerLimit(
                    H,
                    Bd_turn_new
                )
            )

        else:

            # ------------------------------------------------
            # STEP 362
            # CLEAR ISLATCHED
            # ------------------------------------------------

            IsLatched_new = False


            if Bd_new is None:
                raise RuntimeError(
                    "Step 366 requires Bd, "
                    "but Bd is undefined"
                )


            # =================================================
            # STEP 366
            #
            # IS Upper(H) - Bd < Binit(H)
            # AND
            # H <= 0 ?
            # =================================================

            if (
                Upper(H) - Bd_new < Binit(H)
                and H <= 0.0
            ):

                # STEP 368
                # SET ONINITMAG

                OnInitMag_new = True

                # STEP 372
                # B = Binit(H)

                B = Binit(H)

            else:

                # STEP 370
                # B = Upper(H) - Bd

                B = Upper(H) - Bd_new


    # ========================================================
    # RETURN UPDATED STATE
    # ========================================================

    return (
        B,
        OnInitMag_new,
        IsLatched_new,
        Bd_new,
        Bd_turn_new,
        H0_new,
        scale_new,
        offset_new
    )