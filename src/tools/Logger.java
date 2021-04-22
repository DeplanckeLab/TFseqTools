package tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import model.ErrorMessage;

public class Logger 
{
	public static BufferedWriter log = null;
	
	public static void init(File logFilePath)
	{
		if(logFilePath != null)
		{
			try
			{
				Logger.log = new BufferedWriter(new FileWriter(logFilePath));
			}
			catch(IOException ioe)
			{
				new ErrorMessage(ioe.getMessage());
			}
		}
	}
	
	public static void write(String toWrite)
	{
		if(log != null)
		{
			try
			{
				Logger.log.write(toWrite);
			}
			catch(IOException ioe)
			{
				new ErrorMessage(ioe.getMessage());
			}
		}
	}
	
	public static void close()
	{
		if(Logger.log != null)
		{
			try
			{
				Logger.log.close();
			}
			catch(IOException ioe)
			{
				new ErrorMessage(ioe.getMessage());
			}
		}
	}
}
