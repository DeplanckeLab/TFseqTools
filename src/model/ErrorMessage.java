package model;

public class ErrorMessage 
{
	public String displayed_error = "Error";
	
	public ErrorMessage(String errorMessage) 
	{
		if(errorMessage != null) this.displayed_error = errorMessage;
		System.err.println(this.displayed_error);
		System.exit(-1);
	}
}
