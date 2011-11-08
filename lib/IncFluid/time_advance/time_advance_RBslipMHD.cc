#include "../IncFluid.h"


//*********************************************************************************************

void IncFluid::Time_advance(IncVF& W, IncSF& T)
{

	//*********************************************************************************************
	//	EULER
	//

	if (integ_scheme == "EULER") {
		Compute_rhs(W, T);
		Single_time_step(W, T, Tdt, 1, 1, 1);
	}


	//*********************************************************************************************
	//	RK2
	//

	else if (integ_scheme == "RK2") {

		// Allocated once and for all- Vcopy, Wcopy  (static)
		static PlainCVF Vcopy(N);
		static PlainCVF Wcopy(N);
		static PlainCSF Scopy(N);


		Copy_field_to(Vcopy);
		W.Copy_field_to(Wcopy);
		T.Copy_field_to(Scopy);		// Vcopy[i] <- V[i] ; Scopy <- T.F

		Compute_rhs(W, T);

		Single_time_step(W, T, Tdt, 0.5, 0.5, 0.5);  		// Goto the mid point

		Compute_force_TO_rhs(W, T);

		Copy_field_from(Vcopy);
		W.Copy_field_from(Wcopy);
		T.Copy_field_from(Scopy);

		Single_time_step(W, T, Tdt, 1, 0.5, 1);
		// Time-step by Tdt now using the mid-point slopes
	}



	//*********************************************************************************************
	//	RK4
	//

	else if (integ_scheme == "RK4")	{

		static PlainCVF Vcopy(N);
		static PlainCVF Wcopy(N);
		static PlainCSF Scopy(N);

		static PlainCVF tot_Vrhs(N);
		static PlainCVF tot_Wrhs(N);
		static PlainCSF tot_Srhs(N);

		Copy_field_to(Vcopy);
		W.Copy_field_to(Wcopy);
		T.Copy_field_to(Scopy);

		tot_Vrhs.PlainCV_Initialize();
		tot_Wrhs.PlainCV_Initialize();
		tot_Srhs.PlainCS_Initialize();


		// Step 1

		Compute_rhs(W, T);
		Single_time_step(W, T, Tdt, 0.5, 0.5, 0.5);			// u1: Goto the mid point

		Compute_RK4_parts(W, T, tot_Vrhs, tot_Wrhs, tot_Srhs, Tdt, 1.0, Tdt/6);

		// Step 2

		Compute_force_TO_rhs(W, T);

		Copy_field_from(Vcopy);
		W.Copy_field_from(Wcopy);
		T.Copy_field_from(Scopy);

		Single_time_step(W, T, Tdt, 0.5, 0, 0.5);					// u2: Goto the mid pt

		Compute_RK4_parts(W, T, tot_Vrhs, tot_Wrhs, tot_Srhs, Tdt, 0.5, Tdt/3);

		// Step 3

		Compute_force_TO_rhs(W, T);

		Copy_field_from(Vcopy);
		W.Copy_field_from(Wcopy);
		T.Copy_field_from(Scopy);

		Single_time_step(W, Tdt, 1, 0.5, 1);					// u3 : Goto the mid pt

		Compute_RK4_parts(W, T, tot_Vrhs, tot_Wrhs, tot_Srhs, Tdt, 0.5, Tdt/3);

		// Step 4

		Compute_force_TO_rhs(W, T);

		Compute_RK4_parts(W, T, tot_Vrhs, tot_Wrhs, tot_Srhs, Tdt, 0, Tdt/6);

		// Final result

		Copy_field_from(Vcopy);
		Mult_field_exp_ksqr_dt(Tdt, 1.0);

		W.Copy_field_from(Wcopy);
		W.Mult_field_exp_ksqr_dt(Tdt, 1.0);
		T.Copy_field_from(Scopy);
		T.Mult_field_exp_ksqr_dt(Tdt, 1.0);

		*V1 = *V1 + (*tot_Vrhs.V1);
		*V2 = *V2 + (*tot_Vrhs.V2);
		*V3 = *V3 + (*tot_Vrhs.V3);

		*W.V1 = *W.V1 + (*tot_Wrhs.V1);
		*W.V2 = *W.V2 + (*tot_Wrhs.V2);
		*W.V3 = *W.V3 + (*tot_Wrhs.V3);

	}			// Of RK4

	// if free_slip_verticalwall condition is initialized as IC, it should
	// be satisfied at all times, yet we set it again just to make sure everytime step.
	if ((free_slip_verticalwall_switch == 1) && (basis_type == "SCFT"))
		free_slip_verticalwall_field(T);

	//
	// For all schemes
	if (alias_switch == "DEALIAS")		Dealias(W, T);					// Keep V(k) dealiased

}

//*****************************   End of time_advance_RBslipMHD.cc  ***************************
