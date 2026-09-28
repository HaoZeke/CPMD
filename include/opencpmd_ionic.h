/*
 * Ionic surface shared by an external optimizer and one OpenCPMD run.
 *
 * The wavefunction stays allocated for the whole run. An optimizer
 * supplies Cartesian positions in bohr and receives the total energy
 * in hartree and the CPMD ionic forces in hartree/bohr. The force is
 * the negative gradient. rgsaddle's surface callback wants the
 * gradient, so that driver stores minus these forces.
 *
 * This header is the whole published contract. eOn, rgsaddle, and
 * gpr_optim are separate libraries. This tree does not contain them.
 *
 * eOn talks in files. eon_geom is an atom count and then one x y z
 * line per ion in angstrom, in &ATOMS species order. A count of 0, a
 * count that is not the deck, or a missing file ends the run.
 * eon_force is the total energy in hartree and then the forces in
 * hartree/bohr, same ion order.
 *
 * rgsaddle and gpr_optim are loaded at run time.
 * OPENCPMD_RGSADDLE_LIB names the rgsaddle library
 * (librgsaddle.so otherwise). OPENCPMD_GPR_LIB names the gpr_optim
 * library (libgpr_optim.so otherwise). gpr_optim exports
 *
 *     int opencpmd_gpr_run(const opencpmd_ionic_task *task);
 *
 * and calls task->eval on every MPI rank with the same positions.
 * The deck's cell is already fixed. task->cell is nine bohr values,
 * the columns of the CPMD lattice matrix, and may be used or ignored.
 * Atomic numbers are the periodic-table indexes CPMD stored for each
 * species.
 *
 * rgsaddle reads rgsaddle_task or, when that search is the one to
 * run, rgsaddle_mode. Both files are bohr, species order.
 *
 * rgsaddle_task:
 *   n_images
 *   tangent spring projection method
 *   spring_k force_tol max_move memory
 *   n_images blocks of n_atoms lines, x y z
 *
 * The four integers are the rgsaddle enumerators. Nothing here fills
 * in a tolerance the file left out.
 *
 * rgsaddle_mode:
 *   kind method
 *   dr rotation_tol max_rotations krylov_dim force_tol max_move
 *   n_atoms lines of position, x y z
 *   n_atoms lines of the mode, x y z
 *
 * A present rgsaddle_mode file selects the minimum-mode session.
 * Otherwise rgsaddle_task selects the band. The band steps until
 * rgsaddle reports convergence.
 */
#ifndef OPENCPMD_IONIC_H
#define OPENCPMD_IONIC_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define OPENCPMD_IONIC_ABI_MAJOR 1
#define OPENCPMD_IONIC_ABI_MINOR 0

typedef struct opencpmd_ionic_task opencpmd_ionic_task;

/*
 * One geometry. positions and forces each hold 3 * n_atoms doubles.
 * Return 0 on success, negative on failure. Collective: every rank
 * enters with the same positions.
 */
typedef int (*opencpmd_ionic_eval_fn)(int64_t n_atoms, const double *positions,
                                      double *energy, double *forces);

struct opencpmd_ionic_task {
    int32_t major;
    int32_t minor;
    int64_t n_atoms;
    const double *positions;
    const int32_t *atomic_numbers;
    const double *cell;
    opencpmd_ionic_eval_fn eval;
};

/* Exported by the OpenCPMD binary. The task also carries this pointer. */
int opencpmd_ionic_eval(int64_t n_atoms, const double *positions,
                        double *energy, double *forces);

/* Exported by libgpr_optim. Not defined in this tree. */
int opencpmd_gpr_run(const opencpmd_ionic_task *task);

#ifdef __cplusplus
}
#endif

#endif /* OPENCPMD_IONIC_H */
